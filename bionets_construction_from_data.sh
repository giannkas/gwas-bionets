#!/bin/sh
#
# Script: Data Preprocessing and Analysis Pipeline (Nextflow DSL2 version).
# Description: This script performs the following steps:
#   1. Splitting (or subsampling) data and PLINK files generation
#   2. SNP P-value computation
#   3. MAGMA analysis: annotation and gene analysis steps
#   4. Biological network construction (HotNet2, Hierarchical HotNet, dmGWAS,
#      LEAN, SigMod and Heinz)
# Usage: sh scripts/bionets_construction_from_data.sh
#        Fill in every empty variable below before running. To run it as an
#        interactive job on a Slurm cluster (adjust partition, CPUs and memory):
#        srun -p <partition> -c <cpus> --mem <memory_in_MB> --pty bash scripts/bionets_construction_from_data.sh
# Author: Giann Karlo Aguirre-Samboní
# Date: 20/10/2024
#
# General notes
#   - Every `nextflow2 run` call executes a DSL2 pipeline located in `bin_gwas-bionets2/`.
#   - `k`, `plink` and `profile` are shared by several steps; keep them consistent.
#   - Paths may be absolute or relative to the directory the script is run from.
#   - To switch dataset, change the paths/prefixes of each step together: the
#     output folder of one step is the input folder of the next.
 
 
# =============================================================================
# Step 1: Splitting or subsampling data and PLINK files generation
# =============================================================================
# Splits the genotype data into `k` parts and generates PLINK files for each part.
# The step is skipped when k = 1; Steps 2-4 then read the original data directly.
#
# Parameters:
#   bpfiles       Path to the genotype files, given as the common basename without
#                 extension. Matching files are expected: pgen/bed, pvar/bim and
#                 psam/fam (PLINK 2 / PLINK 1 binary formats respectively).
#   k             Number of parts the data is split into. Must be an integer
#                 greater than 1 for the split to take place.
#   base_out_dir  Output location. At least an output filename (basename) must be given.
#   plink         PLINK version used to write the files (1 or 2). MAGMA expects
#                 PLINK 1 formatted files, so use 1 for this step.
#   profile       Nextflow profile that selects the execution settings (executor,
#                 resources, containers), as defined in the pipeline's
#                 nextflow.config (e.g. local or a cluster profile).
 
bpfiles=""
k=5
base_out_dir=""
plink=1
profile="local"
 
if [ "$k" -gt 1 ]; then
  nextflow2 run bin_gwas-bionets2/split_data.nf \
    --bpfile "$bpfiles" \
    --k "$k" \
    --out "$base_out_dir" \
    --plink "$plink" \
    --profile "$profile"
fi
 
 
# =============================================================================
# Step 2: SNP P-value computation
# =============================================================================
# Computes the association between genetic variants (e.g. SNPs) and the phenotype
# of interest (e.g. disease status), for each of the `k` data splits.
#
# Parameters:
#   bpfolder      Path to the parent folder containing the data splits (the
#                 'data_splits' folder created in Step 1). If no splitting took
#                 place, the folder holding the original data; `prefix` then
#                 selects the files to read.
#   base_out_dir  Output location for the SNP P-values. The same folder is reused
#                 by Steps 3 and 4 for MAGMA and network results.
#   prefix        Basename of the genetic data input (the common name of the
#                 PLINK files, without extension).
#   plink         PLINK version used for the association test (1 or 2). It is set
#                 to 2 here, overriding Step 1, to run the association with
#                 covariates in PLINK 2.
#   covars        Path to the covariate file (tab-separated) adjusted for in the
#                 association test, e.g. age and genetic principal components.
#                 Sample IDs must match those in the genotype files.
#   k, profile    As in Step 1.
 
bpfolder=""
base_out_dir=""
prefix=""
plink=1
covars=""
 
nextflow2 run bin_gwas-bionets2/snps_pvalue.nf \
  --bpfolder "$bpfolder" \
  --k "$k" \
  --prefix "$prefix" \
  --plink "$plink" \
  --covars "$covars" \
  --out "$base_out_dir" \
  --profile "$profile"
 
 
# =============================================================================
# Step 3: Basic analysis with MAGMA: annotation and gene analysis
# =============================================================================
# Uses MAGMA v1.10 to (i) annotate SNPs to genes and (ii) aggregate SNP P-values
# into gene-level P-values (gene analysis), accounting for linkage disequilibrium
# through the genotype reference data.
#
# Parameters:
#   window_size   Size, in kilobases, of the window added around each gene when
#                 annotating SNPs. SNPs within this distance of a gene's
#                 boundaries (upstream and downstream) are assigned to that gene.
#   snplocpval    Folder with the SNP location file(s) (SNP ID, chromosome,
#                 base-pair position) and the SNP P-values from Step 2.
#   geneloc       Path to the gene location file (gene ID, chromosome, start,
#                 stop, strand, gene name). Its genome build must match the
#                 genotype data (e.g. NCBI37.3 for build 37, NCBI38 for build 38).
#   geneannot     Folder where the gene annotation file(s) produced by the
#                 MAGMA annotation step are written and read.
#   prefix        Basename of the genetic data input (same as in Step 2).
#   magma         Path to the MAGMA binary.
#   bpfolder, base_out_dir, k, profile   As in Steps 1 and 2.
 
window_size=50
snplocpval=""
geneloc=""
geneannot=""
magma=""
 
nextflow2 run bin_gwas-bionets2/magma_calc.nf \
  --window "$window_size" \
  --k "$k" \
  --snploc_pval "$snplocpval" \
  --gene_loc "$geneloc" \
  --bpfolder "$bpfolder" \
  --gene_annot "$geneannot" \
  --prefix "$prefix" \
  --magma "$magma" \
  --out "$base_out_dir" \
  --profile "$profile"
 
 
# =============================================================================
# Step 4: Construction of biological networks
# =============================================================================
# Builds disease-associated subnetworks (modules) from the gene-level scores of
# Step 3 over a reference interaction network, using six methods:
#   - Heinz   : finds the maximum-scoring connected subnetwork (BioNet package).
#               P-values are modelled with a Beta-uniform mixture (BUM) model and
#               turned into node scores (positive for significant genes, negative
#               otherwise) using the chosen FDR.
#   - SigMod  : finds a strongly interconnected, high-scoring module via an exact
#               graph-cut algorithm. P-values are converted to scores with
#               score = Phi^-1(1 - p), where Phi^-1 is the inverse CDF of the
#               standard normal distribution. Its own parameters (nmax = 300,
#               maxjump = 10, lambda_max = 1) are not exposed by this script,
#               so the SigMod defaults apply.
#   - dmGWAS  : (dmGWAS_3.0) greedy search for dense modules in a PPI network
#               whose nodes are weighted by gene-based P-values, using
#               node weight = Phi^-1(1 - p). Without expression data, only node
#               weights are used (the dmGWAS_1.0/2.X behaviour), controlled by
#               `d` and `r`.
#   - LEAN    : (LEANR package) computes an enrichment P-value for the network
#               neighbourhood ("local subnetwork") of each gene against a
#               background of randomly drawn P-values, then corrects for
#               multiple testing.
#   - HotNet2 : finds subnetworks with more "heat" (gene score) than expected,
#               using an insulated heat-diffusion (influence) matrix of the network.
#   - Hierarchical HotNet : successor of HotNet2; builds a hierarchy of
#               high-scoring subnetworks from a random-walk similarity matrix,
#               assesses significance on permuted data and summarises results
#               across networks and score sets by consensus.
#
# General parameters:
#   net_ref       Reference network from which the methods construct subnetworks:
#                 a tab-separated list of interacting gene pairs. It must be
#                 undirected and simple (no self-loops or duplicated edges), and
#                 gene identifiers must match exactly (case-sensitive) those in
#                 the gene scores file.
#   net_results   Folder where the results of the methods are saved. The gene
#                 scores file (`magma_scores`) is expected to be located here.
#   magma_scores  Filename of the gene-level scores produced in Step 3, after
#                 conversion to the gene identifiers used by the network
#                 (output of magma_calc.nf, named '<prefix>.scores.genes.out_converted.tsv').
#   data_samp     Sampling origin for the `k` splits, as a boolean (1 or 0):
#                 1 when the splits come from the data (genotypes, Steps 1-3);
#                 0 when they come from the gene scores (or P-values) instead.
#   k, profile    As in Step 1.
#
# Method code locations:
#   sigmod_path   Folder with SigMod's internal code (the R functions of the
#                 SigMod_v2 'R' folder; requires igraph >= 1.0.0).
#   hotnet2_path  Folder with HotNet2's internal code (HotNet2.py,
#                 makeHeatFile.py, makeNetworkFiles.py). Set to a real path only
#                 if HotNet2 is run; it is no longer updated by its authors.
#   hhotnet_path  Folder with the cloned Hierarchical HotNet repository (its
#                 'src' scripts). Input files: index-to-gene, edge list and
#                 gene-to-score files, in the same format as for HotNet2.
#
# Method parameters:
#   fdr           (Heinz) False discovery rate used to calibrate node scores from
#                 the BUM model. Lower values give stricter scoring and smaller
#                 modules; higher values give larger modules. BioNet examples use
#                 values from 1e-3 to 1e-14 for small modules, so 0.5 is permissive.
#   d             (dmGWAS) Order of the neighbour genes searched when expanding a
#                 module: 1 = direct neighbours only, 2 = also neighbours of
#                 neighbours. Order 1 or 2 are used in dmGWAS_1.0/2.X, but
#                 dmGWAS_3.0 uses 1 by default.
#   r             (dmGWAS) Minimum score increment required to add a neighbour to
#                 a module (a neighbour is added only if the module score grows
#                 by more than a factor 1 + r). Small r is permissive and may
#                 include unrelated, low-weight nodes; large r is strict and may
#                 miss genes with moderate P-values. Suggested value: 0.1.
#   lfdr_cutoff   (HotNet2) Local FDR parameter used to sparsify the scores:
#                 genes whose P-value is >= this cutoff receive no heat (a score
#                 of zero).
#   beta          (Hierarchical HotNet) Restart probability of the random walk
#                 used to build the similarity matrix (between 0 and 1). Higher
#                 values keep heat closer to its source, favouring smaller,
#                 more local subnetworks.
#   permutations  (Hierarchical HotNet) Number of permuted datasets used to
#                 assess the significance of the hierarchy clusters.
#   repetitions   Number of random repetitions/background samples used by the
#                 method(s) that estimate empirical P-values. LEAN calls this
#                 `n_reps` (default 10000); larger values give more meaningful
#                 empirical P-values at the cost of run time.
 
net_ref=""
net_results=""
magma_scores="gene_scores_converted.tsv"
fdr=0.5
lfdr_cutoff=0.125
data_samp=1
d=1
r=0.1
sigmod_path=""
hotnet2_path=""
hhotnet_path=""
repetitions=100000
permutations=100
beta=0.5
 
nextflow2 run bin_gwas-bionets2/bionets.nf \
  --network "$net_ref" \
  --k "$k" \
  --d_samp "$data_samp" \
  --scores "$magma_scores" \
  --sigmod_path "$sigmod_path" \
  --hotnet2_path "$hotnet2_path" \
  --hhotnet_path "$hhotnet_path" \
  --heinz_fdr "$fdr" \
  --dmgwas_d "$d" \
  --dmgwas_r "$r" \
  --lfdr "$lfdr_cutoff" \
  --reps "$repetitions" \
  --perms "$permutations" \
  --beta "$beta" \
  --out "$net_results" \
  -profile "$profile"
#!/bin/sh
#
# Script: Data Preprocessing and Analysis Pipeline.
# Description: This script performs the following steps:
#   1. Splitting data and PLINK files generation
#   2. SNPs P-value computation
#   3. MAGMA analysis: annotation and gene analysis steps
#   4. Biological network construction
# Usage: ./bionets_construction_from_data.sh
# Author: Giann Karlo Aguirre-Samboní
# Date: 20/10/2024


# Step 1: Splitting data and PLINK files generation
# This step splits the genotype data into `k` parts and generates PLINK files.

# bpfiles: path pattern to pgen/bed, pvar/bim and psam/fam files common name without the extension.
# k: number of times data will be split, it must be an integer gretar than 1.
# base_out_dir: address to save files, at least an output filename must be given.
# plink: version of PLINK to use (1 or 2), however MAGMA expects PLINK file to be formatted in the first version.
# profile: nextflow variable to denote which setting use.
# dsl1: flag to restrict nextflow on domain-specific language 1 (DSL1), pipeline remains to be updated to DSL2.

bpfiles=""
k=5
base_out_dir=""
plink=1
profile=""

if [ $k -gt 1 ]; then
  bin/split_data.nf \
    --bpfile "$bpfiles" \
    --k $k \
    --out "$base_out_dir" \
    --plink $plink \
    -dsl1 \
    -profile "$profile"
fi

# Step 2: SNPs P-value computation
# This step calculates the association between genetic variants (eg. SNPs) and a phenotype of interest (eg. psoriasis).

# bpfolder: path to parent folder where data splits are, named as 'data_splits' or if not splitting process took place, then 
# the folder where data is, prefix parameter will handle to read files.
# prefix: basename of the genetic data input.

bpfolder=""
base_out_dir=""
prefix=""

bin/snps_pvalue.nf \
  --bpfolder "$bpfolder" \
  --k $k \
  --prefix "$prefix" \
  --plink $plink \
  --out "$base_out_dir" \
  -dsl1 \
  -profile "$profile"

# Step 3: Basic analysis with MAGMA: annotation and gene analysis steps
# This step performs an annotation and gene analysis step using MAGMAv1.10 software

# window_size: interspace of SNP influence over a gene of use in MAGMA annotation (size in kilobases).
# snplocpval: folder to where the SNP location files(s) is(are).
# geneloc: path to the gene location file (eg. NCBI37.3.gene.loc).
# geneannot: folder to where the gene annotation files(s) is(are).
# magma: path to magma binary file.

window_size=50
snplocpval=""
geneloc=""
geneannot=""
#prefix=""
prefix=""
magma=""

bin/magma_calc.nf \
  --window "$window_size"\
  --k $k \
  --snploc_pval "$snplocpval" \
  --gene_loc "$geneloc" \
  --bpfolder "$bpfolder" \
  --gene_annot "$geneannot" \
  --prefix "$prefix" \
  --magma "$magma" \
  --plink $plink \
  --out "$base_out_dir" \
  -dsl1 \
  -profile "$profile"

# Step 4: Construction of biological networks using different methods: HotNet2, SigMod and Heinz.

# net_ref: network of reference the methods will use for constructing subnetworks.
# net_results: path to folder to save results out of the methods.
# magma_scores: prefix of the filename where gene scores are stored and it is assumed to be located at net_results folder.
# fdr: false discovery rate for use in Heinz method.
# lfdr_cutoff: local false discovery rate parameter of use for sparse scores of the HotNet2 method (cutoff if P-value >= 0.125).
# data_samp: boolean (1 or 0) to indicate whether the sampling (data splits) comes from the data or the gene scores (or p-values).
# sigmod: path to sigmod files with method's internal code.
# hotnet2: path to hotnet2 files with method's internal code.

net_ref=""
net_results=""
magma_scores=""
fdr=0.5
lfdr_cutoff=0.125
data_samp=1
sigmod_path=""
hotnet2_path=""

bin/bionets.nf \
  --network $net_ref\
  --k $k \
  --d_samp $data_samp \
  --scores $magma_scores \
  --sigmod_path $sigmod_path \
  --hotnet2_path $hotnet2_path \
  --fdr $fdr \
  --lfdr $lfdr_cutoff \
  --out $net_results \
  -dsl1 \
  -profile $profile

# Testing to run MAGMA instead of VEGAS2
# srun -p cbio-cpu -c1 --mem 120000 --pty bash scripts/test_magma.sh


data="/cluster/CBIO/data1/glemoine/ukb_gwas-tools/data/genotypes/basket_4018281/preprocessed/psoriasis/ukb22418_all_chr_b0_v2_filtered_qc-ed"
only_samples="results/psoriasis/magma_scores/snp_association_plink/"

plink --bfile ${data} --assoc --allow-no-sex --out ${only_samples}


# Modifying the plink output to order cols in compatible format for magma
Rscript -e '
  library(dplyr)
  readr::read_table("results/psoriasis/magma_scores/snp_association_plink.tsv") %>% 
  select(SNP, CHR, BP, P) %>%
  readr::write_tsv("results/psoriasis/magma_scores/snp_association_plink_reformat.tsv")
'

# Creating annotation
./scripts/additional_tools/magma/magma \
  --annotate window=50\
  --snp-loc results/psoriasis/magma_scores/snp_association_plink_reformat.tsv \
  --gene-loc /cluster/CBIO/data1/glemoine/ukb_gwas-tools/scripts/additional_tools/magma/genome_data/NCBI37.3/NCBI37.3.gene.loc \
  --out results/psoriasis/magma_scores/magma_output_50kb_window/pso

# Compting score
./scripts/additional_tools/magma/magma \
  --bfile /cluster/CBIO/data1/glemoine/ukb_gwas-tools/data/genotypes/basket_4018281/preprocessed/psoriasis/ukb22418_all_chr_b0_v2_filtered_qc-ed \
  --pval results/psoriasis/magma_scores/snp_association_plink_reformat.tsv N=193608 \
  --gene-annot /cluster/CBIO/data1/gaguirresamboni/ukb_gwas-tools/results/psoriasis/magma_scores/magma_output_50kb_window/pso.genes.annot \
  --out results/psoriasis/magma_scores/magma_output_50kb_window/pso.scores
#--bfile /cluster/CBIO/data1/glemoine/ukb_gwas-tools/scripts/additional_tools/magma/genome_data/g1000_eur/g1000_eur \
# N=193608 from `wc -l data/genotypes/basket_4018281/preprocessed/psoriasis/ukb22418_all_chr_b0_v2_filtered_qc-ed.fam`

# Converting score file gene ids from entrez to hgnc and formating as vegas output for compatibility
Rscript -e '
  library(dplyr)
  magma <- readr::read_table("results/psoriasis/magma_scores/magma_output_50kb_window/pso.scores.genes.out")
  gprofiler2::gconvert(magma$GENE, numeric_ns="ENTREZGENE_ACC") %>% 
    mutate(GENE = as.numeric(input)) %>% 
    select(GENE, name) %>% 
    left_join(magma, by = "GENE") %>%
    rename(Gene = name, Chr = CHR, nSNPs = NSNPS, Start = START, Stop = STOP, Test = ZSTAT, Pvalue = P) %>%
    readr::write_tsv("results/psoriasis/magma_scores/magma_output_50kb_window/pso.scores.genes.out_converted.tsv")
'

# Running subnetwork detection tools
ppi_ref="/cluster/CBIO/data1/glemoine/ukb_gwas-tools/data/network/biogrid_ppi.tsv"
# net_results="results/psoriasis/magma_scores/magma_output_no_window/"
net_results="results/psoriasis/magma_scores/magma_output_50kb_window"
magma_scores="${net_results}/pso.scores.genes.out_converted.tsv"


bin/heinz.nf \
  --vegas $magma_scores \
  --tab2 $ppi_ref \
  --out $net_results/heinz \
  --fdr 0.5 \
  -dsl1 \
  -with-report $net_results/heinz/log \
  -profile cbio_cluster

# bin/lean.nf \
#   --vegas $magma_scores \
#   --tab2 $ppi_ref \
#   --out $net_results/lean \
#   -dsl1 \
#   -with-report $net_results/lean/log \
#   -profile cbio_cluster

bin/sigmod.nf \
  --vegas $magma_scores \
  --tab2 $ppi_ref \
  --out $net_results/sigmod \
  -dsl1 \
  -with-report $net_results/sigmod/log \
  -profile cbio_cluster

bin/hotnet2.nf \
  --scores $magma_scores  \
  --tab2 $ppi_ref \
  --hotnet2_path /gwas-tools/hotnet2 \
  --lfdr_cutoff 0.125 \
  --out $net_results/hotnet2 \
  -dsl1 \
  -with-report $net_results/hotnet2/log \
  -profile cbio_cluster 

# bin/dmgwas.nf \
#   --vegas $magma_scores \
#   --tab2 $ppi_ref \
#   --out $net_results/dmgwas \
#   -dsl1 \
#   -with-report $net_results/dmgwas/log \
#   -profile cbio_cluster


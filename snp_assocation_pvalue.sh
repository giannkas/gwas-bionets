# SNPs P-value computation

bfiles="/cluster/CBIO/data1/gaguirresamboni/ukb_gwas-tools/data/genotypes/basket_4018281/preprocessed/psoriasis"
k=5
base_out_dir="results/psoriasis/magma_scores"
prefix="ukb22418_all_chr_b0_v2_filtered_qc-ed"

bin/snps_pvalue.nf \
  --bfile "$bfiles" \
  --k "$k" \
  --prefix "$prefix" \
  --out "$base_out_dir" \
  -dsl1 \
  -profile cbio_cluster
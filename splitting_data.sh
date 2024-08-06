# splitting data

bfiles="/cluster/CBIO/data1/glemoine/ukb_gwas-tools/data/genotypes/basket_4018281/preprocessed/psoriasis/ukb22418_all_chr_b0_v2_filtered_qc-ed"
k=5
base_out_dir="data/genotypes/basket_4018281/preprocessed/psoriasis"


bin/split_data.nf \
  --bfile "$bfiles" \
  --k "$k" \
  --out "$base_out_dir" \
  -dsl1 \
  -profile cbio_cluster
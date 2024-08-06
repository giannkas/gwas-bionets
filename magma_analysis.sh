# Basic analysis with MAGMA: annotation step and gene analysis step

window_size=50
k=5
snplocpval="/cluster/CBIO/data1/gaguirresamboni/ukb_gwas-tools/results/psoriasis/magma_scores"
geneloc="/cluster/CBIO/data1/glemoine/ukb_gwas-tools/scripts/additional_tools/magma/genome_data/NCBI37.3/NCBI37.3.gene.loc"
base_out_dir="results/psoriasis/magma_scores"
geneannot="/cluster/CBIO/data1/gaguirresamboni/ukb_gwas-tools/results/psoriasis/magma_scores"
bfiles="/cluster/CBIO/data1/gaguirresamboni/ukb_gwas-tools/data/genotypes/basket_4018281/preprocessed/psoriasis"
prefix="ukb22418_all_chr_b0_v2_filtered_qc-ed"
profile="cbio_cluster"

# MAGMA steps: annotation and gene analysis
# ./scripts/additional_tools/magma/magma \
bin/magma_calc.nf \
  --window "$window_size"\
  --k "$k" \
  --snploc_pval "$snplocpval" \
  --gene_loc "$geneloc" \
  --bfile "$bfiles" \
  --prefix "$prefix" \
  --pval "$snplocpval" \
  --gene_annot "$geneannot" \
  --out "$base_out_dir" \
  -dsl1 \
  -profile "$profile"

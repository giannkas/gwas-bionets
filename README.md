# gwas-bionets
_This repository is a based on [hclimente/gwas-tools](https://github.com/hclimente/gwas-tools) which uses the GWAS analysis part to prepare its adaptation for another project at [CBIO](https://cbio.mines-paristech.fr/)_

1. Splitting data if needed. Open the file `splitting_data.sh` to change parameters.

`scripts/splitting_data.sh`

2. SNPs P-value computation using PLINK v1.9. Open the file `snp_association_pvalue.sh` to change parameters. Note that there is a function to reformat PLINK output to order cols in compatible format for MAGMA.

`scripts/snp_association_pvalue.sh`

3. MAGMA analysis to conduct annotation and gene analysis steps. Open the file `magma_analysis.sh` to change parameters.

`scripts/magma_analysis.sh`

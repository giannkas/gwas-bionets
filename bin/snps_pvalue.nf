#!/usr/bin/env nextflow

params.out = '.'
params.k = 5

/////////////////////////////////////
//  SNP P-VALUES COMPUTATION
/////////////////////////////////////

process compute_snps_pvalue {

  publishDir "${params.out}/data_splits/data_split_${task.index}", mode: "copy"

  input:
    val data from params.bfile
    val prx from params.prefix
    val I from 1..params.k
    val K from params.k

  output:
    path "snp_association_plink_split_${I}.assoc" into snppvalues

  script:
    """
      plink --bfile ${data}/data_splits/data_split_${I}/${prx}_split_${I} --assoc --allow-no-sex --out snp_association_plink_split_${I}

    """
}

/////////////////////////////////////
//  SNP P-VALUES REFORMAT
/////////////////////////////////////

process snps_pvalue_reformat {

  publishDir "${params.out}/data_splits/data_split_${task.index}", mode: "copy"

  input:
    val I from 1..params.k
    path assoc from snppvalues

  output:
    path "snp_association_plink_split_${I}.tsv" into new_snppvalues

  script:
    """
      #!/usr/bin/env Rscript

        library(dplyr)
        readr::read_table("${assoc}") %>% 
        select(SNP, CHR, BP, P) %>%
        readr::write_tsv("snp_association_plink_split_${I}.tsv")
    """

}

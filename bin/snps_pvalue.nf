#!/usr/bin/env nextflow

log.info ""
log.info "========================================================"
log.info "|          [gwas-bionets] - snps_pvalue.nf          |"
log.info "========================================================"
log.info ""
log.info "### General association analysis between SNPs and the trait (eg. psoriasis) ###"
log.info ""
log.info ""
log.info ""
log.info "--------------------------------------------------------"
log.info "This program comes with NO WARRANTY"
log.info "It is free software, see LICENSE for details about"
log.info "redistribution and contribution."
log.info "--------------------------------------------------------"
log.info ""

// Help info ##########
params.help = null
if (params.help) {
    log.info ""
    log.info "Usage : snps_pvalue.nf --bpfolder <parent_folder> --k <knumber> --prefix <my_prefix> --plink <pversion> --out <filename>"
    log.info ""
    log.info "  --bpfolder    path to parent folder where genotypes has been split and where 'data_splits' folder resides so it"
    log.info "                can search iteratively the genetic data input (pgen/bed, pvar/bim and psam/fam files)"
    log.info "                for each data split (1, 2, ..., k)."
    log.info "  --k           number of times data was split, it must be in accordance to the previous k for splitting the data."
    log.info "  --prefix      prefix or basename for the genetic data input (eg. my_prefix.bed, my_prefix.bim, my_prefix.fam)."
    log.info "  --plink       version of PLINK to use (1 or 2), note that pgen, pvar and psam is for PLINK2, whereas"
    log.info "                bed, bim and fam corresponds to PLINK."
    log.info "  --out         path/to/filename where output files will be saved."
    log.info ""
    log.info ""
    log.info ""
    log.info "Example : "
    log.info "snps_pvalue.nf \\"
    log.info "  --bpfolder           path/to/parent_folder \\"
    log.info "  --k                 5 \\"
    log.info "  --prefix            my_prefix \\"
    log.info "  --plink             1 \\"
    log.info "  --out               path/to/my_output \\"
    log.info "  -dsl1 \\"
    log.info "  -profile            my_cluster \\"
    log.info ""

    exit 0
}


params.out = '.'
params.k = 1
params.plink = 1

/////////////////////////////////////
//  SNP P-VALUES COMPUTATION
/////////////////////////////////////

process compute_snps_pvalue {

  publishDir { params.out + (params.k > 1 ? '/data_splits/data_split_' + task.index : '') }, overwrite: true, mode: "copy"

  input:
    val data from params.bpfolder
    val prx from params.prefix
    val I from 1..params.k
    val K from params.k
    val V from params.plink


  output:
    file { "snp_association_plink" + (K > 1 ? "_split_" + I : '') + (V == 1 ? ".assoc" : ".PHENO1.glm.logistic.hybrid") } into snppvalues

  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""
    def dir_splits = K > 1 ? "/data_splits/data_split_${I}" : ""

    """
    if [ $V -eq 1 ]; then
      plink --bfile ${data}${dir_splits}/${prx}${split_suffix} --assoc --allow-no-sex --out snp_association_plink${split_suffix}
    else
      plink2 --pfile ${data}${dir_splits}/${prx}${split_suffix} --glm allow-no-covars --out snp_association_plink${split_suffix}
    fi

    """
}

/////////////////////////////////////
//  SNP P-VALUES REFORMAT
//
//  Modifying PLINK output to order 
//  cols in compatible format for magma
/////////////////////////////////////

process snps_pvalue_reformat {

  publishDir { params.out + (params.k > 1 ? '/data_splits/data_split_' + task.index : '') }, overwrite: true, mode: "copy"


  input:
    val I from 1..params.k
    path assoc from snppvalues
    val K from params.k
    val V from params.plink

  output:
    path { "snp_association_plink" + (K > 1 ? "_split_" + I : '') + ".tsv" } into new_snppvalues

  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""

    """
      #!/usr/bin/env Rscript

        library(dplyr)
        if (${V} == 1) {
          readr::read_table("${assoc}") %>% 
          select(SNP, CHR, BP, P) %>%
          readr::write_tsv("snp_association_plink${split_suffix}.tsv")
        } else {
          readr::read_table("${assoc}") %>% 
          select(ID, `#CHROM`, POS, P) %>%
          rename(SNP = ID, CHR = `#CHROM`, BP = POS) %>%
          readr::write_tsv("snp_association_plink${split_suffix}.tsv")
        }
        
    """

}

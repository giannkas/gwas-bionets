#!/usr/bin/env nextflow

/////////////////////////////////////
//  SNP P-VALUES COMPUTATION
/////////////////////////////////////

process compute_snps_pvalue {

  publishDir { params.out + ((params.k as Integer) > 1 ? '/stable_consensus/data_splits/data_split_' + I : '') }, overwrite: true, mode: "copy"

  input:
    path data
    val prefix
    val K
    val vplink
    path covars
    each I


  output:
    path { "${prefix}_snp_association_plink" + (K > 1 ? "_split_" + I : '') + (vplink == 1 ? ".assoc" : ".PHENO1.glm.logistic.hybrid") }, emit: snppvalues

  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""
    def dir_splits = K > 1 ? "/data_splits/data_split_${I}" : ""
    // println "K=${K} I=${I}"
    //plink2 --pfile ${data}${dir_splits}/${prefix}${split_suffix} --glm allow-no-covars hide-covar --out ${prefix}_snp_association_plink${split_suffix}
    """
    if [ $vplink -eq 1 ]; then
      plink --bfile ${data}${dir_splits}/${prefix}${split_suffix} --assoc --allow-no-sex --out ${prefix}_snp_association_plink${split_suffix}
    else
      plink2 --bfile ${data}${dir_splits}/${prefix}${split_suffix} --glm --covar ${covars} --covar-name sex age PC1 PC2 --covar-variance-standardize --out ${prefix}_snp_association_plink${split_suffix}
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

  publishDir { params.out + ((params.k as Integer) > 1 ? '/stable_consensus/data_splits/data_split_' + task.index : '') }, overwrite: true, mode: "copy"


  input:
    path snppvalues
    val prefix
    val K
    val vplink

  output:
    path { "${prefix}_snp_association_plink" + (K > 1 ? "_split_" + task.index : '') + ".tsv" }, emit: new_snppvalues

  script:
    def split_suffix = K > 1 ? "_split_${task.index}" : ""

    """
      #!/usr/bin/env Rscript

        library(dplyr)
        if ($vplink == 1) {
          readr::read_table("$snppvalues") %>% 
          select(SNP, CHR, BP, P) %>%
          readr::write_tsv("${prefix}_snp_association_plink${split_suffix}.tsv")
        } else {
          readr::read_table("$snppvalues") %>% 
          filter(TEST == "ADD") %>%
          select(ID, `#CHROM`, POS, P) %>%
          rename(SNP = ID, CHR = `#CHROM`, BP = POS) %>%
          readr::write_tsv("${prefix}_snp_association_plink${split_suffix}.tsv")
        }
        
    """

}

params.out = '.'
params.k = 1
params.plink = 1
params.covars = ""
params.bpfolder = ""
params.prefix = ""

// Help info ##########
params.help = null

workflow {
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

  def K = params.k as Integer
  def vplink = params.plink as Integer

  // Step 1: Compute SNP p-values
  def snppvalues = compute_snps_pvalue(
    params.bpfolder, params.prefix, K, vplink, params.covars, 1..K)

  // Step 2: Reformat SNP p-values
  def new_snppvalues = snps_pvalue_reformat(
      snppvalues, params.prefix, K, vplink)

}
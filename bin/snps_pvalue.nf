#!/usr/bin/env nextflow

// =============================================================================
// snps_pvalue.nf  -  Step 2 of gwas-bionets: SNP-level association P-values
// =============================================================================
// Tests the association between each SNP and the trait (e.g. disease status) in
// every data split, and writes the results in the format MAGMA expects.
// The phenotype is read from the sample file (.fam, or .psam for PLINK 2).
//
// Workflow
//   1. compute_snps_pvalue   PLINK association test.
//        --plink 1: allelic case/control test (plink --assoc), output '.assoc'
//        --plink 2: logistic regression (plink2 --glm), output '.PHENO1.glm.logistic.hybrid'
//   2. snps_pvalue_reformat  R script that keeps SNP, CHR, BP and P in a tab-separated
//                            file (for PLINK 2 only the additive-effect rows, TEST == ADD).
//
// Parameters (command line, as --name value)
//   --bpfolder  Parent folder of the genotype data. With k > 1 the files of split I are read
//               from <bpfolder>/data_splits/data_split_<I>/<prefix>_split_<I>; with k = 1 from
//               <bpfolder>/<prefix> (the --out of split_data.nf).
//   --prefix    Basename of the genetic data input (e.g. my_prefix.bed/.bim/.fam).
//   --k         Number of splits; it must match the --k used in split_data.nf (default 1).
//   --plink     PLINK version, 1 (default) or 2.
//   --covars    Covariates file for the PLINK 2 regression (see NOTE below).
//   --out       Output folder (default '.').
//   --help      Print the usage and exit.
//
// Output layout
//   k > 1 : <out>/stable_consensus/data_splits/data_split_<I>/<prefix>_snp_association_plink_split_<I>.{tsv,...}
//   k = 1 : <out>/<prefix>_snp_association_plink.{tsv,...}
//   The raw PLINK result is published next to the reformatted '.tsv' (columns SNP, CHR, BP, P).
//
// Requirements: 'plink' / 'plink2' in the PATH; R with the packages dplyr and readr.

/////////////////////////////////////
//  SNP P-VALUES COMPUTATION
/////////////////////////////////////

// compute_snps_pvalue: one task per split I (1..K).
//   data     --bpfolder: parent folder of the genotype data
//   prefix   --prefix
//   K, vplink   number of splits (--k) and PLINK version (--plink)
//   covars   --covars: covariates file
//   I        split index (1..K)
//   emits    snppvalues: '<prefix>_snp_association_plink[_split_<I>]' + '.assoc' (PLINK 1)
//            or '.PHENO1.glm.logistic.hybrid' (PLINK 2)
//
// NOTE on covariates (current behaviour of the shell block below):
//   - '\$covars' is a SHELL variable, which is never defined (Nextflow only substitutes
//     '${covars}'), so the test [ "\$covars" = "" ] is always true: with --plink 2 the
//     'allow-no-covars' branch always runs and --covars is NOT applied.
//   - The last branch (covariates) would call plink2 with --bfile and fixed covariate names
//     (sex, age, PC1, PC2); those columns must then exist in the covariates file.

process compute_snps_pvalue {

  publishDir { params.out + ((params.k as Integer) > 1 ? '/stable_consensus/data_splits/data_split_' + I : '') }, overwrite: true, mode: "copy"

  input:
    path data
    val prefix
    val K
    val vplink
    val covars          // an absolute path, or "" for none
    each I

  output:
    path { "${prefix}_snp_association_plink" + (K > 1 ? "_split_" + I : '') + (vplink == 1 ? ".assoc" : ".PHENO1.glm.logistic.hybrid") }, emit: snppvalues

  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""
    def dir_splits   = K > 1 ? "/data_splits/data_split_${I}" : ""
    def infile       = "${data}${dir_splits}/${prefix}${split_suffix}"
    def outname      = "${prefix}_snp_association_plink${split_suffix}"

    // decided by Nextflow, not by the shell
    def has_covars = covars.toString().trim() != ''
    def glm_opts   = has_covars ?
        "--glm hide-covar --covar ${covars} --covar-name sex age PC1 PC2 --covar-variance-standardize" :
        "--glm allow-no-covars hide-covar"

    """
    if [ $vplink -eq 1 ]; then
      plink --bfile ${infile} --assoc --allow-no-sex --out ${outname}
    else
      # split files may be PLINK 1 (bed/bim/fam) or PLINK 2 (pgen/pvar/psam)
      if [ -f ${infile}.pgen ]; then fmt="--pfile"; else fmt="--bfile"; fi
      plink2 \$fmt ${infile} ${glm_opts} --out ${outname}
    fi
    """
}

/////////////////////////////////////
//  SNP P-VALUES REFORMAT
//
//  Modifying PLINK output to order 
//  cols in compatible format for magma
/////////////////////////////////////

// snps_pvalue_reformat: one task per split. Converts the PLINK result into a
// tab-separated file with columns SNP, CHR, BP, P, as needed by MAGMA (--snp-loc, --pval).
//   PLINK 1: columns SNP, CHR, BP, P are kept as they are.
//   PLINK 2: keeps TEST == 'ADD' rows and renames ID, #CHROM, POS to SNP, CHR, BP.
//   emits    new_snppvalues: '<prefix>_snp_association_plink[_split_<n>].tsv'

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

params.out = '.'             // --out: output folder
params.k = 1                 // --k: number of splits (same as in split_data.nf)
params.plink = 1             // --plink: PLINK version (1: --assoc, 2: --glm)
params.covars = ""           // --covars: covariates file (PLINK 2; see NOTE in compute_snps_pvalue)
params.bpfolder = ""         // --bpfolder: parent folder of the genotype data (see header)
params.prefix = ""           // --prefix: basename of the genetic data input

// Help info ##########
params.help = null            // --help: print the usage and exit

workflow {
  log.info ""
  log.info "========================================================"
  log.info "|          [gwas-bionets] - snps_pvalue.nf          |"
  log.info "========================================================"
  log.info ""
  log.info "### General association analysis between SNPs and the trait (eg. disease status) ###"
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
    log.info "Usage : snps_pvalue.nf --bpfolder <parent_folder> --k <knumber> --prefix <my_prefix> --plink <pversion> [--covars <file>] --out <filename>"
    log.info ""
    log.info "  --bpfolder    path to parent folder where genotypes has been split and where 'data_splits' folder resides so it"
    log.info "                can search iteratively the genetic data input (pgen/bed, pvar/bim and psam/fam files)"
    log.info "                for each data split (1, 2, ..., k)."
    log.info "  --k           number of times data was split, it must be in accordance to the previous k for splitting the data."
    log.info "  --prefix      prefix or basename for the genetic data input (eg. my_prefix.bed, my_prefix.bim, my_prefix.fam)."
    log.info "  --plink       version of PLINK to use (1 or 2, default 1), note that pgen, pvar and psam is for PLINK2, whereas"
    log.info "                bed, bim and fam corresponds to PLINK. 1: allelic association test (--assoc); 2: logistic"
    log.info "                regression (--glm)."
    log.info "  --covars      covariates file used by the PLINK 2 regression (default none). See the NOTE in the process"
    log.info "                compute_snps_pvalue: it is currently not applied."
    log.info "  --out         path/to/folder where output files will be saved. With k > 1 results go to"
    log.info "                <out>/stable_consensus/data_splits/data_split_<I>/."
    log.info "  --help        show this message and exit."
    log.info ""
    log.info "Example :"
    log.info "snps_pvalue.nf \\"
    log.info "  --bpfolder          path/to/parent_folder \\"
    log.info "  --k                 5 \\"
    log.info "  --prefix            my_prefix \\"
    log.info "  --plink             1 \\"
    log.info "  --out               path/to/my_output \\"
    log.info "  -profile            my_cluster \\"

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
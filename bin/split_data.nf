#!/usr/bin/env nextflow

/////////////////////////////////////
//  SPLIT PREPARATION
/////////////////////////////////////

process make_splits {

  publishDir { params.out + ((params.k as Integer) > 1 ? '/data_splits/data_split_' + I : '') }, overwrite: true, mode: "copy"

  input:
    path genbed
    path samfam
    val K
    each I
    val vplink

  output:
    path {"${genbed.baseName}" + (K > 1 ? "_split_" + I : '') + (vplink == 1 ? ".fam" : ".psam") }, emit: splits

  script:    
    template 'genotypes/split_samples.sh'
  

}

///////////////////////////////////////////////////////////////////////////////////
//  GENOTYPES, SAMPLES AND VARIANT INFORMATION GENERATION FROM A SAMPLE FILE
///////////////////////////////////////////////////////////////////////////////////

process genbed_varbim_samfam_creation {

  publishDir { params.out + ((params.k as Integer) > 1 ? '/data_splits/data_split_' + task.index : '') }, mode: "copy"

  input:
    val data
    path splits
    path genbed
    val K
    val vplink

  output:
    path {"${genbed.baseName}" + (K > 1 ? "_split_" + task.index : '') + (vplink == 1 ? ".bim" : ".pvar") }, emit: varbims
    path {"${genbed.baseName}" + (K > 1 ? "_split_" + task.index : '') + (vplink == 1 ? ".bed" : ".pgen") }, emit: genbeds
    path {"${genbed.baseName}" + (K > 1 ? "_split_" + task.index : '') + (vplink == 1 ? ".fam" : ".psam") }, emit: samfams

  script:
    def split_suffix = K > 1 ? "_split_${task.index}" : ""

    """
    if [ $vplink -eq 1 ]; then
      plink --bfile ${data} --keep ${splits} --make-bed --out ${genbed.baseName}${split_suffix}
    else
      plink2 --pfile ${data} --keep ${splits} --make-pgen --out ${genbed.baseName}${split_suffix}
    fi
    """
}

// Help info ##########
params.help = null

// default values for the parameters
params.out = '.'
params.k = 1
params.plink = 1
params.bpfile = null

workflow {
  log.info ""
  log.info "========================================================"
  log.info "|          [gwas-bionets] - split_data.nf          |"
  log.info "========================================================"
  log.info ""
  log.info "### Splitting samples data and generating new corresponding PLINK files ###"
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
    log.info "Usage : split_data.nf --bpfile <filepattern> --k <knumber> --out <filename>"
    log.info ""
    log.info "  --bpfile      path to pgen/bed, pvar/bim and psam/fam files and only the file prefix needs to be specified."
    log.info "  --k           number of times data will be split, it must be a non-negative integer."
    log.info "  --out         path/to/filename where output files will be saved."
    log.info ""
    log.info ""
    log.info ""
    log.info "Example : "
    log.info "split_data.nf \\"
    log.info "  --bpfile             path/to/bpfiles \\"
    log.info "  --k                 5 \\"
    log.info "  --out               path/to/myoutput \\"
    log.info "  -dsl1 \\"
    log.info "  -profile            my_cluster \\"
    log.info ""

    exit 0
  }

  // Input Validation
  if (!params.bpfile) {
    log.error "Error: --bpfile parameter must be specified."
    exit 1
  }

  // File assignments based on PLINK version
  def genbed
  def varbim
  def samfam
  def vplink = params.plink as Integer

  if (vplink == 1) {
    genbed = file("${params.bpfile}.bed")
    varbim = file("${genbed.baseName}.bim")
    samfam = file("${params.bpfile}.fam")
  } else {
    genbed = file("${params.bpfile}.pgen")
    varbim = file("${genbed.baseName}.pvar")
    samfam = file("${params.bpfile}.psam")
  }

  //println "The value of samfam is: ${samfam}"
  def K = params.k as Integer

  // Step 1: Prepare sample splits
  def splits = make_splits(
    genbed, 
    samfam, 
    K, 
    1..K, 
    vplink
  )

  // Step 2: Generate genotype, variant, and sample files
  def outputs = genbed_varbim_samfam_creation(params.bpfile, splits, genbed, K, vplink)

}
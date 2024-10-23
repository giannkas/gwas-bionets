#!/usr/bin/env nextflow

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

// Help info ##########
params.help = null
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

// default values of the parameters.
params.out = '.'
params.k = 1
params.plink = 1

// gwas
genbed = ""
varbim = ""
samfam = ""

// we assign the variables differently depeding on the PLINK version.
if (params.plink == 1) {
  genbed = file("${params.bpfile}.bed")
  varbim = file("${genbed.baseName}.bim")
  samfam = file("${params.bpfile}.fam")
 } else {
  genbed = file("${params.bpfile}.pgen")
  varbim = file("${genbed.baseName}.pvar")
  samfam = file("${params.bpfile}.psam")
 }

/////////////////////////////////////
//  SPLIT PREPARATION
/////////////////////////////////////
process make_splits {

  //publishDir "${params.out}/data_splits/data_split_${task.index}", overwrite: true, mode: "copy"
  publishDir { params.out + (params.k > 1 ? '/data_splits/data_split_' + task.index : '') }, overwrite: true, mode: "copy"


  input:
    file SAMPLES from samfam
    val I from 1..params.k
    val K from params.k
    val V from params.plink

  output:
    path {"${genbed.baseName}" + (K > 1 ? "_split_" + I : '') + (V == 1 ? ".fam" : ".psam") } into splits

  script:
  // the script is saved in the 'template' folder whithin 'genotypes' subfolder.
  // ex. split_data_folder/templates/genotypes/split_samples.sh
  template 'genotypes/split_samples.sh'
  

}


///////////////////////////////////////////////////////////////////////////////////
//  GENOTYPES, SAMPLES AND VARIANT INFORMATION GENERATION FROM A SAMPLE FILE
///////////////////////////////////////////////////////////////////////////////////

process genbed_varbim_samfam_creation {

  publishDir { params.out + (params.k > 1 ? '/data_splits/data_split_' + task.index : '') }, mode: "copy"


  input:
    val data from params.bpfile
    val I from 1..params.k
    val K from params.k
    path insamfams from splits
    val V from params.plink

  output:
    path {"${genbed.baseName}" + (K > 1 ? "_split_" + I : '') + (params.plink == 1 ? ".bim" : ".pvar") } into varbims
    path {"${genbed.baseName}" + (K > 1 ? "_split_" + I : '') + (params.plink == 1 ? ".bed" : ".pgen") } into genbeds
    path {"${genbed.baseName}" + (K > 1 ? "_split_" + I : '') + (params.plink == 1 ? ".fam" : ".psam") } into samfams

  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""

    """
    if [ $V -eq 1 ]; then
      plink --bfile ${data} --keep ${insamfams} --make-bed --out ${genbed.baseName}${split_suffix}
    else
      plink2 --pfile ${data} --keep ${insamfams} --make-pgen --out ${genbed.baseName}${split_suffix}
    fi
    """
}


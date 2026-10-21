#!/usr/bin/env nextflow

// =============================================================================
// split_data.nf  -  Step 1 of gwas-bionets: split samples and write PLINK files
// =============================================================================
// Splits the samples of a PLINK genotype dataset into K parts and writes a new
// PLINK fileset (genotypes, variants and samples) for each part.
//
// Workflow
//   1. make_splits                    Builds, for each split I = 1..K, the list of
//                                     samples that belong to it. The splitting logic
//                                     lives in the template 'genotypes/split_samples.sh'
//                                     (folder 'templates' next to this script).
//   2. genbed_varbim_samfam_creation  Runs PLINK with --keep to extract those samples
//                                     and write the genotype, variant and sample files.
//
// Parameters (command line, as --name value)
//   --bpfile  (required) Path to the input fileset WITHOUT extension. Expected files:
//             .bed/.bim/.fam (PLINK 1) or .pgen/.pvar/.psam (PLINK 2).
//   --k       Number of splits (default 1). Use an integer > 1 to split the data;
//             with 1 the output carries no split suffix and no 'data_splits' folders.
//   --out     Output folder (default '.').
//   --plink   PLINK version, 1 (default) or 2. MAGMA (Step 3) needs PLINK 1 files.
//   --help    Print the usage and exit.
//
// Output layout (<basename> = input fileset name without folder or extension)
//   k > 1 : <out>/data_splits/data_split_<I>/<basename>_split_<I>.{bed,bim,fam}
//           (or .{pgen,pvar,psam} with --plink 2)
//   k = 1 : <out>/<basename>.{bed,bim,fam}  (or .{pgen,pvar,psam})
//
// Requirements: 'plink' (PLINK 1.9) or 'plink2' in the PATH, according to --plink.

/////////////////////////////////////
//  SPLIT PREPARATION
/////////////////////////////////////

// make_splits: one task per split I (1..K).
//   genbed   genotype file of the input fileset (.bed/.pgen); only its basename is used
//   samfam   sample file of the input fileset (.fam/.psam)
//   K        total number of splits (--k)
//   I        index of the split handled by the task (1..K)
//   vplink   PLINK version (--plink)
//   emits    'splits': sample file '<basename>[_split_<I>].fam|.psam' listing the
//            samples assigned to split I (used by --keep in the next process)
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

// genbed_varbim_samfam_creation: one task per split, extracts the split's samples.
//   data     path of the input fileset without extension (--bpfile)
//   splits   sample file of the split, from make_splits
//   genbed, K, vplink   as above
//   PLINK 1: plink  --bfile <data> --keep <splits> --make-bed  --out <basename>[_split_<n>]
//   PLINK 2: plink2 --pfile <data> --keep <splits> --make-pgen --out <basename>[_split_<n>]
//   emits    varbims (.bim|.pvar), genbeds (.bed|.pgen), samfams (.fam|.psam)
//   The split number <n> is Nextflow's task.index (1..K).
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
params.help = null            // --help: print the usage and exit

// default values for the parameters
params.out = '.'             // --out: output folder
params.k = 1                 // --k: number of splits (> 1 to split the samples)
params.plink = 1             // --plink: PLINK version (1: bed/bim/fam, 2: pgen/pvar/psam)
params.bpfile = null         // --bpfile: input fileset path without extension (required)

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
    log.info "Usage : split_data.nf --bpfile <filepattern> --k <knumber> --out <filename> [--plink <pversion>]"
    log.info ""
    log.info "  --bpfile      (required) path/to/fileset, without extension. The pgen/pvar/psam (PLINK 2) or bed/bim/fam (PLINK 1)"
    log.info "                files must share this prefix."
    log.info "  --k           number of parts the samples are split into, an integer greater than 1 (default 1: no splitting,"
    log.info "                files are written without split suffix or 'data_splits' folders)."
    log.info "  --out         path/to/folder where output files will be saved (default '.'). With k > 1 each split is written"
    log.info "                to <out>/data_splits/data_split_<I>/ (I = 1..k)."
    log.info "  --plink       version of PLINK to use (1 or 2, default 1). Note that pgen, pvar and psam is for PLINK2, whereas"
    log.info "                bed, bim and fam corresponds to PLINK. MAGMA (next steps) expects PLINK 1 files."
    log.info "  --help        show this message and exit."
    log.info ""
    log.info "Example :"
    log.info "split_data.nf \\"
    log.info "  --bpfile            path/to/bpfiles \\"
    log.info "  --k                 5 \\"
    log.info "  --out               path/to/myoutput \\"
    log.info "  --plink             1 \\"
    log.info "  -profile            my_cluster \\"

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
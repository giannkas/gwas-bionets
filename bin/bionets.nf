#!/usr/bin/env nextflow

//include { dmgwas_nf as dmgwas } from './dmgwas'
//include { heinz_nf as heinz } from './heinz'
//include { hotnet2_nf as hotnet2 } from './hotnet2.nf'
//include { sigmod_call as sigmod } from './sigmod'


// General parameters
params.out = '.'
params.k = 1
params.d_samp = 1
//params.profile = 'local'

// SigMod parameters
params.sigmod_lambdamax = 1
params.sigmod_nmax = 300
params.sigmod_maxjump = 10

// Heinz parameters
params.heinz_fdr = 0.2

// dmGWAS parameters
params.dmgwas_r = 0.1
params.dmgwas_d = 2

// Help info ##########
params.help = null

/////////////////////////////////////
//  HEINZ METHOD
/////////////////////////////////////

process heinz_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val fdr
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/heinz.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --tab2 ${network} \
        --out ${out}${dir_splits}/heinz \
        --fdr ${fdr} \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}

/////////////////////////////////////
//  DMGWAS METHOD
/////////////////////////////////////

process dmgwas_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val d
      val r
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/dmgwas.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --tab2 ${network} \
        --out ${out}${dir_splits}/dmgwas \
        --d ${d} \
        --r ${r} \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}


/////////////////////////////////////
//  SIGMOD METHOD
/////////////////////////////////////

process sigmod_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val sigmod_path
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/sigmod.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --network ${network} \
        --sigmod_path ${sigmod_path} \
        --out ${out}${dir_splits}/sigmod \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}

/////////////////////////////////////
//  HOTNET2 METHOD
/////////////////////////////////////

process hotnet2_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val lfdr
      val hotnet2_path
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/hotnet2.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --tab2 ${network} \
        --lfdr_cutoff ${lfdr} \
        --hotnet2_path ${hotnet2_path} \
        --out ${out}${dir_splits}/hotnet2 \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}

/////////////////////////////////////
//  HIERARCHICAL HOTNET METHOD
/////////////////////////////////////

process hhotnet_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val conn
      val perms
      val beta
      val hhotnet_path
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/hhotnet.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --network ${network} \
        --connectivity ${conn} \
        --permutations ${perms} \
        --beta ${beta} \
        --hhotnet_path ${hhotnet_path} \
        --out ${out}${dir_splits}/hhotnet \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}

/////////////////////////////////////
//  LEAN METHOD
/////////////////////////////////////

process lean_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val reps
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/lean.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --network ${network} \
        --reps ${reps} \
        --out ${out}${dir_splits}/lean \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}


// Workflow
workflow {
  log.info ""
  log.info "========================================================"
  log.info "|          [gwas-bionets] - bionets.nf          |"
  log.info "========================================================"
  log.info ""
  log.info "### Script to call different biological network methods ###"
  log.info ""
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
    log.info "Usage : bionets.nf --network <network_file> --k <knumber> --d_samp <data_sampling> --scores <scores_file> \\"
    log.info "                      --fdr <false_discovery_rate> --lfdr <lfdr_cutoff>  --out <filename>\\"
    log.info ""
    log.info ""
    log.info "  --network     network of reference for the methods to operate and build a network according to gene scores (p-values)."
    log.info "  --k           number of times data was split, it must be in accordance to the previous k for splitting the data."
    log.info "  --d_samp      1 or 0 (default 1) to specify where it was sampled from. For instance, you sample the data in different splits or"
    log.info "                one can divide the gene scores file in different chunks and use each chunk as your gene scores (or p-values) file."
    log.info "  --scores      file for the formatted gene scores, produced by the 'snps_pvalue_reformat' process of the magma_calc.nf file."
    log.info "  --fdr         false discovery rate parameter to control the resultant subnetwork size in the Heinz method."
    log.info "  --lfdr        local fdr of use in HotNet2 method to avoid false positives in multiple testing analyses."
    log.info "                "
    log.info "                "
    log.info ""
    log.info ""
    log.info ""
    log.info "Example : "
    log.info "bionets.nf \\"
    log.info "  --network           path/to/my_network \\"
    log.info "  --k                 5 \\"
    log.info "  --d_samp            1 \\"
    log.info "  --scores            path/to/my_scores \\"
    log.info "  --bpfolder          path/to/genotypes_folder \\"
    log.info "  --fdr               0.5 \\"
    log.info "  --lfdr              0.125 \\"
    log.info "  --out               path/to/my_output \\"
    log.info "  -dsl2 \\"
    log.info "  -profile            my_cluster \\"
    log.info ""

    exit 0
  }

  // Define inputs
  def scores = file(params.scores)
  def network = file(params.network)
  def sigmod_path = params.sigmod_path
  def hotnet2_path = params.hotnet2_path
  def hhotnet_path = params.hhotnet_path
  def K = params.k as Integer
  def d_samp = params.d_samp as Integer
  def fdr = params.heinz_fdr as Double
  def lfdr = params.lfdr as Double
  def conn =  params.conn as Integer
  def perms = params.perms as Integer
  def beta = params.beta as Double
  def D = params.dmgwas_d as Integer
  def R = params.dmgwas_r as Double
  def lambdamax = params.sigmod_lambdamax as Integer
  def nmax = params.sigmod_nmax as Integer
  def maxjump = params.sigmod_maxjump as Integer
  def reps = params.reps as Integer

  // println "The value of scores is: ${params.scores}"
  // println "The value of out is: ${params.out}"
  // println "The value of scores baseName is: ${scores.baseName}"

  // Run processes inline
  heinz_call(scores, network, K, d_samp, params.out, fdr, 1..K)
  dmgwas_call(scores, network, K, d_samp, params.out, D, R, 1..K)
  sigmod_call(scores, network, K, d_samp, params.out, sigmod_path, 1..K)
  hotnet2_call(scores, network, K, d_samp, params.out, lfdr, hotnet2_path, 1..K)
  lean_call(scores, network, K, d_samp, params.out, reps, 1..K)
  // hhotnet_call(scores, network, K, d_samp, params.out, conn, perms, beta, hhotnet_path, 1..K)
  
}


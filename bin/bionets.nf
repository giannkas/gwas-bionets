#!/usr/bin/env nextflow

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

// Help info ##########
params.help = null
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
    log.info "  -dsl1 \\"
    log.info "  -profile            my_cluster \\"
    log.info ""

    exit 0
}


params.out = '.'
params.k = 1
params.d_samp = 1
params.profile = params.profile ?: 'local'


/////////////////////////////////////
//  HEINZ METHOD
/////////////////////////////////////

process heinz_call {

  input:
    val magma from params.scores
    val net from params.network
    val fdr from params.fdr
    val I from 1..params.k
    val K from params.k
    val samp from params.d_samp

  script:
    def split_suffix = (K > 1 && samp != 0) ? "_split_${I}" : (K > 1 && samp == 0) ? "_chunk_${I}" : ""
    def dir_splits = (K > 1 && samp != 0) ? "/data_splits/data_split_${I}" : (K > 1 && samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""
    def ith = ( K > 1 ) ? I : 0

    """
      ${baseDir}/heinz.nf \
        --scores ${params.out}${dir_splits}/${magma}${split_suffix}.tsv \
        --tab2 ${net} \
        --out ${params.out}${dir_splits}/heinz \
        --fdr ${fdr} \
        --i ${ith} \
        --d_samp ${samp} \
        -dsl1 \
        -with-report ${params.out}${dir_splits}/heinz/log \
        -profile ${params.profile}
    """
}

/////////////////////////////////////
//  SIGMOD METHOD
/////////////////////////////////////

process sigmod_call {

  input:
    val magma from params.scores
    val net from params.network
    val I from 1..params.k
    val K from params.k
    val samp from params.d_samp
    val sigmod from params.sigmod_path

    script:
      def split_suffix = (K > 1 && samp != 0) ? "_split_${I}" : (K > 1 && samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (K > 1 && samp != 0) ? "/data_splits/data_split_${I}" : (K > 1 && samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""
      def ith = ( K > 1 ) ? I : 0

    """
      ${baseDir}/sigmod.nf \
        --scores ${params.out}${dir_splits}/${magma}${split_suffix}.tsv \
        --tab2 ${net} \
        --sigmod_path ${sigmod}
        --out ${params.out}${dir_splits}/sigmod \
        --i ${ith} \
        --d_samp ${samp} \
        -dsl1 \
        -with-report ${params.out}${dir_splits}/sigmod/log \
        -profile ${params.profile}
    """

}

/////////////////////////////////////
//  HOTNET2 METHOD
/////////////////////////////////////

process hotnet2_call {

  input:
    val magma from params.scores
    val net from params.network
    val I from 1..params.k
    val lfdr from params.lfdr
    val K from params.k
    val samp from params.d_samp
    val hotnet2 from params.hotnet2_path

  script:
    def split_suffix = (K > 1 && samp != 0) ? "_split_${I}" : (K > 1 && samp == 0) ? "_chunk_${I}" : ""
    def dir_splits = (K > 1 && samp != 0) ? "/data_splits/data_split_${I}" : (K > 1 && samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""
    def ith = ( K > 1 ) ? I : 0

    """
      ${baseDir}/hotnet2.nf \
        --scores ${params.out}${dir_splits}/${magma}${split_suffix}.tsv  \
        --tab2 ${net} \
        --lfdr_cutoff ${lfdr} \
        --hotnet2_path ${hotnet2}
        --out ${params.out}${dir_splits}/hotnet2 \
        --i ${ith} \
        --d_samp ${samp} \
        -dsl1 \
        -with-report ${params.out}${dir_splits}/hotnet2/log \
        -profile ${params.profile}
    """

}

#!/usr/bin/env nextflow

process sigmod {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path SCORES
    path NETWORK
    val SIGMOD_PATH
    val LAMBDA_MAX
    val NMAX
    val MAXJUMP
    val SPLIT

  output:
    path "selected_genes${SPLIT}.sigmod.txt", emit: genes_sigmod

  script:
    template 'discovery/run_sigmod.R'
}

params.out = '.'
params.lambdamax = 1
params.nmax = 500
params.maxjump = 30
params.sigmod_path = null
params.i = 0
params.d_samp = 1

workflow {

  def img_sigmod = "/gwas-bionets/sigmod"
  def sigmod_path = (params.sigmod_path != null && params.sigmod_path != "/default/path") ? params.sigmod_path : img_sigmod

  // Split suffix
  def split = (params.i > 0 && params.d_samp != 0) ? "_split_${params.i}" : 
            (params.i > 0 && params.d_samp == 0) ? "_chunk_${params.i}" : ""

  println "The value of sigmod_path is: ${sigmod_path}"

  // input files
  def scores = file(params.scores)
  def network = file(params.network)

  // run sigmod process
  sigmod(scores, network, sigmod_path, params.lambdamax, params.nmax, params.maxjump, split)
}

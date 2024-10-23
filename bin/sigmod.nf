#!/usr/bin/env nextflow

params.out = '.'
params.lambdamax = 1
params.nmax = 300
params.maxjump = 10
params.sigmod_path = null
docker_sigmod = '/gwas-bionets/sigmod'
params.i = 0
params.d_samp = 1


split = params.i > 0 && params.d_samp != 0 ? "_split_${params.i}" : params.i > 0 && params.d_samp == 0 ? "_chunk_${params.i}" : ""


// conditional SigMod input to handle Docker or Dockerless execution
if (params.sigmod_path != null) {
    SIGMOD_PATH = file(params.sigmod_path) 
} else {
    SIGMOD_PATH = docker_sigmod
}

// annotation
MAGMA_OUT = file(params.scores)
TAB2 = file(params.tab2)

process sigmod {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file MAGMA_OUT
        file TAB2
        if (params.sigmod_path != null) {
           file SIGMOD_PATH
        } else {
            val SIGMOD_PATH
        }
        val LAMBDAMAX from params.lambdamax
        val NMAX from params.nmax
        val MAXJUMP from params.maxjump
        val SPLIT from split

    output:
        file "selected_genes${SPLIT}.sigmod.txt" into genes_sigmod

    script:
    template 'discovery/run_sigmod.R'
}

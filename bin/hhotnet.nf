#!/usr/bin/env nextflow

process compile_fortran {
    publishDir "${HHOTNET}/src", mode: 'copy'

    input:
      path HHOTNET
    
    output:
      path "fortran_module*.so"
    
    script:
      """    
      f2py -c ${HHOTNET}/src/fortran_module.f95 -m fortran_module
      """
}

process data_formatting {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path NETWORK
    path SCORES

  output:
    path 'index_gene.tsv', emit: index_gene
    path 'gene_score.tsv', emit: gene_score
    path 'edge_list.tsv', emit: edge_list

  script:
    template 'io/scrs_net_format_hhotnet.R'

}

process similarity_matrix {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path HHOTNET
    path NETWORK
    val BETA

  output:
    path "similarity_matrix.h5", emit: sim_mtx
    path "beta.txt", emit: beta_file

  script:
    """
      python ${HHOTNET}/src/construct_similarity_matrix.py \
        -i   ${NETWORK} \
        -o   similarity_matrix.h5 \
        -bof beta.txt
    """
}

process permuting_network {

  input:
    tuple val(I), path(NETWORK), path(HHOTNET)

  output:
    path "edge_list_${I}.tsv", emit: edge_list_$I

  script:
    """
      python ${HHOTNET}/src/permute_network.py \
        -i ${NETWORK} \
        -s "${I}" \
        -c \
        -o edge_list_${I}.tsv
    """

}

process find_permutation_bins {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path NETWORK
    path SCORES
    path IDXGENE
    path HHOTNET

  output:
    path "score_bins.tsv", emit: score_bins

  script:
    """
      python ${HHOTNET}/src/find_permutation_bins.py \
        -gsf ${SCORES} \
        -igf ${IDXGENE} \
        -elf ${NETWORK} \
        -ms 1000 \
        -o score_bins.tsv
    """

}

process permuting_scores {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    tuple val(PERMS), path(BINS), path(SCORES), path(HHOTNET)

  output:
    tuple val(PERMS), path("scores_${PERMS}.tsv")

  script:
    """
      python ${HHOTNET}/src/permute_scores.py \
        -i ${SCORES} \
        -bf ${BINS} \
        -s "${PERMS}" \
        -o scores_${PERMS}.tsv
    """

}

process construct_hierarchies {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    tuple val(PERMS), path(SIMMATRIX), path(IDXGENE), path(SCORES), path(HHOTNET)

  output:
    tuple val(PERMS), path("hierarchy_edge_list_${PERMS}.tsv"), path("hierarchy_index_gene_${PERMS}.tsv")

  script:
    """
      python ${HHOTNET}/src/construct_hierarchy.py \
        -smf ${SIMMATRIX} \
        -igf ${IDXGENE} \
        -gsf ${SCORES} \
        -helf hierarchy_edge_list_${PERMS}.tsv \
        -higf hierarchy_index_gene_${PERMS}.tsv
    """

}

process processing_hierarchies {

  publishDir params.out, overwrite: true, mode: 'copy'
  debug true
  
  input:
    path OBS_HEDGE_LIST
    path OBS_HINDEX_GENE
    path PERM_HEDGE_LISTS
    path PERM_HINDEX_GENES
    val NETWORKNAME
    val SCORESNAME
    path HHOTNET
    val SPLIT

  output:
    path "clusters_${NETWORKNAME}_${SCORESNAME}${SPLIT}.hhotnet.tsv", emit: clusters
    path "sizes_${NETWORKNAME}_${SCORESNAME}${SPLIT}.hhotnet.pdf", emit: sizes

  script:
    """
      python ${HHOTNET}/src/process_hierarchies.py \
        -oelf ${OBS_HEDGE_LIST} \
        -oigf ${OBS_HINDEX_GENE} \
        -pelf ${PERM_HEDGE_LISTS} \
        -pigf ${PERM_HINDEX_GENES} \
        -lsb 10 \
        -cf clusters_${NETWORKNAME}_${SCORESNAME}${SPLIT}.hhotnet.tsv \
        -pl ${NETWORKNAME} gene_score \
        -pf sizes_${NETWORKNAME}_${SCORESNAME}${SPLIT}.hhotnet.pdf
    """

}

process performing_consensus {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path CLUSTERS
    path IDXGENE
    path NETWORK
    val NETWORKNAME
    val SCORESNAME
    path HHOTNET
    val SPLIT

  output:
    path "selected_genes${SPLIT}.hhotnet.tsv"
    path "consensus_edges${SPLIT}.hhotnet.tsv"

  script:
    """
      python ${HHOTNET}/src/perform_consensus.py \
        -cf ${CLUSTERS} \
        -igf ${IDXGENE} \
        -elf ${NETWORK} \
        -n ${NETWORKNAME} \
        -s ${SCORESNAME} \
        -t 1 \
        -cnf selected_genes${SPLIT}.hhotnet.tsv \
        -cef consensus_edges${SPLIT}.hhotnet.tsv
    """
}


// Parameters
params.out = "."
params.hhotnet_path = null
params.i = 0
params.d_samp = 1
params.connectivity = 1
params.beta = 0.5
params.permutations = 5

workflow {
  // input files
  def network = file(params.network)
  def scores = file(params.scores)

  def split = (params.i > 0 && params.d_samp != 0) ? "_split_${params.i}" : 
            (params.i > 0 && params.d_samp == 0) ? "_chunk_${params.i}" : ""

  def img_hhotnet = "/gwas-bionets2/hhotnet"
  def HHOTNET = (params.hhotnet_path != null && params.hhotnet_path != "/default/path") ? params.hhotnet_path : img_hhotnet

  // def network_permutations = 100
  // def heat_permutations = 1000
  def network_permutations = params.permutations
  def beta = params.beta
  def conn = params.connectivity
  def network_name = network.getBaseName()
  def scores_name = scores.getBaseName()

  // processes
  data_formatting(network, scores)
  // folders(network, scores)
  similarity_matrix(HHOTNET, data_formatting.out.edge_list, beta)
  net_nperms = Channel.of(1..4)

  net_perms = data_formatting.out.edge_list
    .combine(net_nperms)
    .map { edge_list, nperm -> [nperm, edge_list, HHOTNET] }

  // permuted_nets = permuting_network(net_perms)


  // permutation bins

  find_permutation_bins(
    data_formatting.out.edge_list,
    data_formatting.out.gene_score,
    data_formatting.out.index_gene,
    HHOTNET)

  // permutation scores

  chl1_nperms = Channel.of(1..network_permutations)

  scores_perms = find_permutation_bins.out.score_bins
    .combine(data_formatting.out.gene_score)
    .combine(chl1_nperms)
    .map { score_bins, gene_score, nperm -> [nperm, score_bins, gene_score, HHOTNET] }

  permuted_scores = permuting_scores(scores_perms)

  // hierarchies construction
  obs_hierarchy_input = Channel.of(0)
    .combine(similarity_matrix.out.sim_mtx)
    .combine(data_formatting.out.index_gene)
    .combine(data_formatting.out.gene_score)
    .map { perm_id, sim_mtx, idx_gene, gene_scr ->
        tuple(perm_id, sim_mtx, idx_gene, gene_scr, HHOTNET)
    }

  perm_hierarchy_input = permuted_scores
    .combine(similarity_matrix.out.sim_mtx)
    .combine(data_formatting.out.index_gene)
    .map { perm_id, gene_scr, sim_mtx, idx_gene ->
        tuple(perm_id, 
              sim_mtx, 
              idx_gene, 
              gene_scr, 
              HHOTNET)
    }

  // combine observed and permuted inputs
  hierarchies_cns = obs_hierarchy_input.mix(perm_hierarchy_input)

  cnstd_hierarchies = construct_hierarchies(hierarchies_cns)

  // hierarchies processing
  // separate observed (perm_id=0) from permuted hierarchies
  obs_hierarchy = cnstd_hierarchies
    .filter { perm_id, hedge_list, hindex_gene -> perm_id == 0 }
    .map { perm_id, hedge_list, hindex_gene -> 
        tuple(hedge_list, hindex_gene) 
  }
  
  perm_hierarchies = cnstd_hierarchies
    .filter { perm_id, hedge_list, hindex_gene -> perm_id > 0 }
    .map { perm_id, hedge_list, hindex_gene -> 
        tuple(hedge_list, hindex_gene) 
  }

  // collect permuted hierarchies
  perm_hedge_lists = perm_hierarchies
    .map { hedge_list, hindex_gene -> hedge_list }
    .collect()
  
  perm_hindex_genes = perm_hierarchies
    .map { hedge_list, hindex_gene -> hindex_gene }
    .collect()

  // process hierarchies
  processing_hierarchies(hierarchies_procs)*/
  processing_hierarchies(
    obs_hierarchy.map { it[0] },  // obs_hedge_list
    obs_hierarchy.map { it[1] },  // obs_hindex_gene
    perm_hedge_lists,
    perm_hindex_genes,
    network_name,
    scores_name,
    HHOTNET,
    split
  )

  // perform consensus
  performing_consensus(processing_hierarchies.out.clusters,
    data_formatting.out.index_gene,
    data_formatting.out.edge_list,
    network_name,
    scores_name,
    HHOTNET,
    split)

}
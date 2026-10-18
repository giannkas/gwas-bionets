#!/usr/bin/env nextflow

process make_network {

  input:
    path TAB2
    path MAGMA

  output:
    path 'node_index.tsv', emit: node_index
    path 'edge_list.tsv', emit: edge_list

  script:
    template 'io/tab2_2hotnet.R'

}

process sparse_scores {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    file MAGMA
    val LFDR_CUTOFF
    val SPLIT

  output:
    path "scored_genes${SPLIT}.sparse.txt", emit: sparse_scores
    path "lfdr_plot${SPLIT}.pdf"

  """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(twilight)
    library(cowplot)

    theme_set(theme_cowplot())

    scores <- read_tsv('${MAGMA}') %>%
      select(Gene, Pvalue) %>%
      filter(Gene != "None") %>%
      mutate(Gene = toupper(Gene)) %>%
      group_by(Gene) %>%
      filter(Pvalue == min(Pvalue)) %>%
      ungroup() %>%
      distinct(Gene, .keep_all = TRUE)

    lfdr <- twilight(scores\$Pvalue, B=1000)
    lfdr <- tibble(Gene = scores\$Gene[as.numeric(rownames(lfdr\$result))],
                magma_p = scores\$Pvalue[as.numeric(rownames(lfdr\$result))],
                lfdr = lfdr\$result\$fdr)

    ggplot(lfdr, aes(x = magma_p, y = 1 - lfdr)) +
     geom_line() +
     geom_vline(xintercept = ${LFDR_CUTOFF}, color = 'red') +
     labs(x = 'P-value', y = '1 - lFDR')
    ggsave('lfdr_plot${SPLIT}.pdf', width=7, height=6)

    lfdr %>%
     mutate(Pvalue = ifelse(magma_p < ${LFDR_CUTOFF}, magma_p, 1)) %>%
     write_tsv('scored_genes${SPLIT}.sparse.txt')
  """

}

process magma2hotnet {

  input:
    path MAGMA

  output:
    path 'scores.ht', emit: scores

  script:
    template 'io/magma2hotnet.R'

}


process make_h5_network {

  input:
    path HOTNET2
    path NODE_INDEX
    path EDGE_LIST
    val BETA
    val NETWORK_PERMUTATIONS

  output:
    path "ppin_ppr_${BETA}.h5", emit: h5
    path 'permuted', emit: permutations

  script:

    """
      python2 ${HOTNET2}/makeNetworkFiles.py \
        --edgelist_file ${EDGE_LIST} \
        --gene_index_file ${NODE_INDEX} \
        --network_name ppin \
        --prefix ppin \
        --beta ${BETA} \
        --cores -1 \
        --num_permutations ${NETWORK_PERMUTATIONS} \
        --output_dir .
    """

 }

process make_heat_data {

  input:
    path HOTNET2
    path SCORES

  output:
    path 'heat.json', emit: heat

  script:
    """
    python2 ${HOTNET2}/makeHeatFile.py \
      scores \
      --heat_file ${SCORES} \
      --output_file heat.json \
      --name gwas
    """
}

process hotnet2 {

  input:
    val HOTNET2
    path HEAT
    path NETWORK
    path PERMUTATIONS
    val NETWORK_PERMUTATIONS
    val HEAT_PERMUTATIONS
    val BETA

  output:
    path 'consensus/subnetworks.tsv', emit: subnetworks

  script:
    """
    python2 ${HOTNET2}/HotNet2.py \
      --network_files ${NETWORK} \
      --permuted_network_path ${PERMUTATIONS}/ppin_ppr_${BETA}_##NUM##.h5 \
      --heat_files ${HEAT} \
      --network_permutations ${NETWORK_PERMUTATIONS} \
      --heat_permutations ${HEAT_PERMUTATIONS} \
      --num_cores -1 \
      --output_directory .
    """

}

process process_output {

  publishDir params.out, overwrite: true, mode: 'copy'

  input:
    path SUBNETWORKS
    val SPLIT

  output:
    path "selected_genes${SPLIT}.hotnet2.tsv", emit: genes_hotnet2

  script:
    """
    #!/usr/bin/env Rscript

     library(tidyverse)

     read_tsv('${SUBNETWORKS}', col_types = 'cc', comment = '#', col_names = F) %>%
         select(X1) %>%
         mutate(cluster = 1:n()) %>%
         separate_rows(X1, sep = ' ') %>%
         rename(gene = X1) %>%
         write_tsv('selected_genes${SPLIT}.hotnet2.tsv')
    """

}

// Parameters
params.out = "."
params.hotnet2_path = null
params.i = 0
params.d_samp = 1
params.lfdr_cutoff = 0.05

workflow {
  // input files
  def tab2 = file(params.tab2)
  def magma = file(params.scores)

  def split = (params.i > 0 && params.d_samp != 0) ? "_split_${params.i}" : 
            (params.i > 0 && params.d_samp == 0) ? "_chunk_${params.i}" : ""

  def img_hotnet2 = "/gwas-bionets/hotnet2"
  def HOTNET2 = (params.hotnet2_path != null && params.hotnet2_path != "/default/path") ? params.hotnet2_path : img_hotnet2

  def network_permutations = 100
  def heat_permutations = 1000
  // def network_permutations = 10
  // def heat_permutations = 20
  def beta = 0.4

  // processes
  make_network(tab2, magma)
  sparse_scores(magma, params.lfdr_cutoff, split)
  magma2hotnet(sparse_scores.out.sparse_scores)

  make_h5_network(HOTNET2, make_network.out.node_index, make_network.out.edge_list, beta, network_permutations)
  make_heat_data(HOTNET2, magma2hotnet.out.scores)

  hotnet2(HOTNET2, make_heat_data.out.heat, make_h5_network.out.h5, 
    make_h5_network.out.permutations, network_permutations, 
    heat_permutations, beta)
  process_output(hotnet2.out.subnetworks, split)
}
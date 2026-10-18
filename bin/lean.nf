#!/usr/bin/env nextflow

process lean {

  publishDir params.out, overwrite: true, mode: "copy"

    input:
      path SCORES
      path EDGELIST
      val REPS
      val SPLIT

    output:
      path "selected_genes${SPLIT}.lean.txt"

    script:

    """
      #!/usr/bin/env Rscript

      library(igraph)
      library(LEANR)
      library(dplyr)
      library(tidyverse)

      # gene scores
      gene_scores <- read_tsv("${SCORES}") %>%
        select(Gene, Pvalue) %>%
        filter(Gene != "None") %>%
        mutate(Gene = toupper(Gene)) %>%
        group_by(Gene) %>%
        filter(Pvalue == min(Pvalue)) %>%
        ungroup() %>%
        distinct(Gene, .keep_all = TRUE)

      # read edgelist
      network <- read_tsv("${EDGELIST}") %>%
        rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
        mutate(gene1 = toupper(gene1), gene2 = toupper(gene2)) %>%
        filter(gene1 %in% gene_scores\$Gene & gene2 %in% gene_scores\$Gene) %>%
        filter(gene1 != gene2) %>% # this line could remove singleton scored genes in the network interacting with themselves.
        filter(!is.na(gene1) & !is.na(gene2)) %>%
        select(gene1, gene2) %>%
        graph_from_data_frame(directed = FALSE)

      # gene_scores <- filter(gene_scores, Gene %in% names(V(network))) 
      # force MAGMA to only include genes in the network
      # this will remove singleton genes, i.e. genes with no connections but which may have a strong signal.

      scores <- gene_scores[['Pvalue']]
      names(scores) <- gene_scores[['Gene']]

      results <- run.lean(scores, network, n_reps = ${REPS}, add.scored.genes = TRUE, ncores = $task.cpus)

      significant_genes <- tibble(Gene = rownames(results\$restab[results\$restab[,'PLEAN']<=0.05,]))

      if (nrow(significant_genes) == 0) {
        warning("No significant genes found. Creating empty output file.")
        tibble(gene = character(0)) %>%
          write_tsv('selected_genes${SPLIT}.lean.txt')
      } else {
        nh_genes <- significant_genes %>%
          mutate(neighbors = map(Gene, ~ {
            if (.x %in% names(results\$nhs)) {
              results\$nhs[[.x]]
            } else {
              character(0)
            }
          })) %>%
          unnest(neighbors)


        nh_genes_unique <- nh_genes %>% 
          rename(gene = neighbors) %>%
          distinct(gene) %>%
          write_tsv('selected_genes${SPLIT}.lean.txt')

      }

    """

}

params.out = '.'
params.reps = 100000
params.network = ''
params.scores = ''

workflow {

  def net = file(params.network)
  def scores = file(params.scores)
  def reps = params.reps
  def split = params.i > 0 && params.d_samp != 0 ? "_split_${params.i}" : params.i > 0 && params.d_samp == 0 ? "_chunk_${params.i}" : ""

  lean(scores, net, reps, split)
}
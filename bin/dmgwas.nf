#!/usr/bin/env nextflow

process dmgwas {

  publishDir params.out, overwrite: true, mode: "copy"

  input:
    path MAGMA
    path EDGELIST
    val D
    val R
    val SPLIT

  output:
    path "selected_genes${SPLIT}.dmgwas.txt"

  """
  #!/usr/bin/env Rscript

  library(dmGWAS)
  library(readr)
  library(dplyr)

  scores <- read_tsv('${MAGMA}') %>% 
    select(Gene, Pvalue) %>%
    mutate(Pvalue = ifelse(Pvalue == 1, 0.99999, Pvalue)) %>%
    filter(Gene != "None") %>%
    mutate(Gene = toupper(Gene)) %>%
    group_by(Gene) %>%
    filter(Pvalue == min(Pvalue)) %>%
    ungroup() %>%
    distinct(Gene, .keep_all = TRUE) %>%
    as.data.frame

  net <- read_tsv('${EDGELIST}') %>%
    rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
    mutate(gene1 = toupper(gene1), gene2 = toupper(gene2)) %>%
    filter(gene1 %in% scores\$Gene & gene2 %in% scores\$Gene) %>%
    filter(gene1 != gene2) %>% # this line could remove singleton scored genes in the network interacting with themselves.
    filter(!is.na(gene1) & !is.na(gene2)) %>%
    select(gene1, gene2) %>%
    as.data.frame

  modules <- dms(net, scores, expr1 = NULL, expr2 = NULL, r = ${R}, d = ${D})
  top <- simpleChoose(modules)

  tibble(gene = names(V(top\$subnetwork))) %>%
      write_tsv("selected_genes${SPLIT}.dmgwas.txt")
  """

}

params.out = '.'
params.r = 0.1
params.d = 2
params.tab2 = ''
params.scores = ''

workflow {

  def tab2 = file(params.tab2)
  def scores = file(params.scores)
  def split = params.i > 0 && params.d_samp != 0 ? "_split_${params.i}" : params.i > 0 && params.d_samp == 0 ? "_chunk_${params.i}" : ""

  dmgwas(scores, tab2, params.d, params.r, split)
}
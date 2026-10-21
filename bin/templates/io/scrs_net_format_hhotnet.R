#!/usr/bin/env Rscript

# NETWORK
# SCORES
# index_gene.tsv, edge_list.tsv, gene_score.tsv

library(tidyverse)
library(igraph)

scores <- read_tsv('${SCORES}') %>%
  select(Gene, Pvalue) %>%
  filter(Gene != "None") %>%
  mutate(Gene = toupper(Gene)) %>%
  group_by(Gene) %>%
  filter(Pvalue == min(Pvalue)) %>%
  ungroup() %>%
  distinct(Gene, .keep_all = TRUE)

network <- read_tsv("${NETWORK}") %>%
  rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
  mutate(gene1 = toupper(gene1), gene2 = toupper(gene2)) %>%
  filter(gene1 %in% scores\$Gene & gene2 %in% scores\$Gene) %>%
  filter(gene1 != gene2) %>% # this line could remove singleton scored genes in the network interacting with themselves.
  filter(!is.na(gene1) & !is.na(gene2)) %>%
  select(gene1, gene2) %>%
  graph_from_data_frame(directed = FALSE)

# Select only MAGMA genes that are also in the PPI - using igraph data structure, Gene-to-score file
gen_scr <- filter(scores, Gene %in% names(V(network)))

# Index-to-gene file
idx_gen <- tibble(id = V(network), names = names(V(network)))

# Edge list file
edg_lst <- as_edgelist(network, names = FALSE) %>% as.data.frame

write_tsv(idx_gen, "index_gene.tsv", col_names = FALSE)
write_tsv(gen_scr, "gene_score.tsv", col_names = FALSE)
write_tsv(edg_lst, "edge_list.tsv", col_names = FALSE)
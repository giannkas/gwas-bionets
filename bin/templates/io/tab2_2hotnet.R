#!/usr/bin/env Rscript
# TAB2
# node_index.tsv, edge_list.tsv
library(tidyverse)
library(igraph)

scores <- read_tsv('${MAGMA}') %>%
  select(Gene, Pvalue) %>%
  filter(Gene != "None") %>%
  mutate(Gene = toupper(Gene)) %>%
  group_by(Gene) %>%
  filter(Pvalue == min(Pvalue)) %>%
  ungroup() %>%
  distinct(Gene, .keep_all = TRUE)

net <- read_tsv("${TAB2}") %>%
  rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
  mutate(gene1 = toupper(gene1), gene2 = toupper(gene2)) %>%
  filter(gene1 %in% scores\$Gene & gene2 %in% scores\$Gene) %>%
  filter(gene1 != gene2) %>% # this line could remove singleton scored genes in the network interacting with themselves.
  filter(!is.na(gene1) & !is.na(gene2)) %>%
  select(gene1, gene2) %>%
  graph_from_data_frame(directed = FALSE)

as_edgelist(net, names = FALSE) %>% as.data.frame %>% write_tsv('edge_list.tsv', col_names = FALSE)
tibble(id = V(net), names = names(V(net))) %>% write_tsv('node_index.tsv', col_names = FALSE)

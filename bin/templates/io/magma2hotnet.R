#!/usr/bin/env Rscript
# MAGMA
# scores.ht

library(tidyverse)

read_tsv('${MAGMA}') %>%
  mutate(score = -log10(Pvalue)) %>%
  select(Gene, score) %>%
  filter(Gene != "None") %>%
  mutate(Gene = toupper(Gene)) %>%
  group_by(Gene) %>%
  filter(score == min(score)) %>%
  ungroup() %>%
  distinct(Gene, .keep_all = TRUE) %>%
  write_tsv('scores.ht', col_names = FALSE)

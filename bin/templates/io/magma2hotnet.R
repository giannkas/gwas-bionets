#!/usr/bin/env Rscript
# MAGMA
# scores.ht

library(tidyverse)

read_tsv('${MAGMA}') %>%
  mutate(score = -log10(Pvalue)) %>%
  select(Gene, score) %>%
  write_tsv('scores.ht', col_names = FALSE)

#!/usr/bin/env Rscript

library(tidyverse)

# load sigmod
require(igraph)
scripts <- list.files('${SIGMOD_PATH}/R', pattern='*.R\$', full.names=TRUE, ignore.case=TRUE)
sapply(scripts, source, .GlobalEnv)
cat("${SIGMOD_PATH}/R: ", "\n", sep="")

# read MAGMA output
scores <- read_tsv('${SCORES}') %>% 
    rename(gene = Gene, p = Pvalue) %>%
    select(gene, p) %>%
    filter(gene != "None") %>%
    mutate(gene = toupper(gene)) %>%
    group_by(gene) %>%
    filter(p == min(p)) %>%
    ungroup() %>%
    distinct(gene, .keep_all = TRUE) %>%
    as.data.frame

# read network
net <- read_tsv("${NETWORK}") %>%
    rename(Interactor_A = `Official Symbol Interactor A`, Interactor_B = `Official Symbol Interactor B`) %>%
    mutate(Interactor_A = toupper(Interactor_A), Interactor_B = toupper(Interactor_B)) %>%
    filter(Interactor_A %in% scores\$gene & Interactor_B %in% scores\$gene) %>%
    filter(Interactor_A != Interactor_B) %>% # this line could remove singleton scored genes in the network interacting with themselves.
    filter(!is.na(Interactor_A) & !is.na(Interactor_B)) %>%
    select(Interactor_A, Interactor_B) %>%
    as.data.frame

# check weight_index = NULL
scored_net <- construct_scored_net(net, interaction_indices = c(1,2), gene_ps = scores)
res_info <- SigMod_bisection(net = scored_net, lambda_max = $LAMBDA_MAX, nmax = $NMAX, maxjump = $MAXJUMP)

save(scored_net, res_info, file = 'sigmod.RData')

data.frame(gene = names(V(res_info\$opt_module[[1]]))) %>%
    write_tsv('selected_genes${SPLIT}.sigmod.txt')

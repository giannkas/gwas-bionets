# =============================================================================
# stable_consensus_net.R  -  Stable consensus network across data splits
# =============================================================================
# Like consensus_net.R, but uses the results of the network methods on every data
# split (Step 4 of gwas-bionets, bionets.nf) to keep only genes that are selected
# repeatedly: a gene is part of the stable consensus when it appears in at least 7 of
# the 25 solutions (5 methods x 5 splits).
# Author: Giann Karlo Aguirre-Samboní
#
# The general procedure of this script is as follows:
#   1. Loads a reference PPI network.
#   2. For every split j = 1..nsplits and every method (Heinz, HotNet2, SigMod, dmGWAS,
#      LEAN), reads the selected genes. It keeps (a) the union of genes per method over the
#      splits, and (b) the number of solutions (method x split) in which each gene appears.
#   3. Keeps the genes found in at least 7 solutions (gene_set7; threshold hard-coded) and
#      annotates them with cGWAS / MAGMABH membership (genes annotated to classical GWAS
#      hits; genes significant in MAGMA after Benjamini-Hochberg correction).
#   4. Writes the gene table and the layout, and saves a figure in which each node is a pie
#      with one slice per method that selected the gene in ANY split; the label fill shows
#      whether the gene is also in MAGMABH and/or cGWAS.
#
# Inputs (edit the paths below; '<prefix>' is the basename used in the pipeline)
#   PPI file       Tab-separated file whose first two columns are the official gene
#                  symbols of each interaction (BioGRID tab2: 'Official Symbol Interactor A/B').
#   results_gwas   Folder with '<prefix>.genes.annot_gene_ids_names.tsv': table with a
#                  'Gene' column listing the genes annotated to the classical GWAS hits.
#   results_44249  Folder '<out>/stable_consensus/data_splits/' of Step 4 (k splits, d_samp = 1),
#                  holding for j = 1..nsplits
#                    data_split_<j>/heinz/selected_genes_split_<j>.heinz.txt
#                    data_split_<j>/hotnet2/selected_genes_split_<j>.hotnet2.tsv  (two columns)
#                    data_split_<j>/sigmod/selected_genes_split_<j>.sigmod.txt
#                    data_split_<j>/dmgwas/selected_genes_split_<j>.dmgwas.txt
#                    data_split_<j>/lean/selected_genes_split_<j>.lean.txt
#                  (each with a 'gene' column).
#   results        Step 4 output folder (--out of bionets.nf), with
#                    <prefix>.scores.genes.out_converted_pruned_bhcorrection.tsv
#                  MAGMA gene results after BH correction and pruning (deduplication): 12 columns, the
#                  third one being 'Gene' (character).
#   nsplits, nmethods  Number of data splits (must match --k) and of network methods (5).
#
# Outputs (<stable> = the 'stable_consensus' folder, i.e. results_44249 + '../')
#   <stable>/stable_consensus_gset7.tsv               Genes of the stable consensus with
#                                                     magmabh_cgwas, cgwas_genes and magmabh flags.
#   <stable>/connections_stable_consensus_gset7.tsv   Layout (x, y, name, xend, yend).
#   <stable>/figures/stable_consensus_methods_gset7.pdf   Figure (the 'figures' folder must exist).
#
# Notes
#   - Only the threshold of 7 is applied (gene_set7), but the user can change that by increasing the number
#     in gene_lists_subsets.
#   - The per-source plots (heinz_plt, hotnet2_plt, ...) are built for individual inspection if needed.
 
# Requirements: R packages igraph, tidygraph, magrittr, ggrepel, tidyverse, cowplot, network,
#   sna, ggnetwork, scatterpie and ggnewscale. ggupset and RCy3 are loaded but not used here.
#   To install tidyverse on Ubuntu you may first need:
#     apt install libcurl4-openssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev


library(igraph)
library(tidygraph)
library(magrittr)
library(ggrepel)
library(tidyverse)
library(cowplot)
suppressMessages(library(network))
suppressMessages(library(sna))
library(ggnetwork)
library(scatterpie)

theme_transparent <- theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
                           plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                           panel.grid.major = element_blank(), # get rid of major grid
                           panel.grid.minor = element_blank(), # get rid of minor grid
                           legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
                           legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
)

theme_set(theme_cowplot() + theme_transparent)

results_gwas <- 'results/psoriasis/magma/magma_output_50kb_window_snp-array_significant_snps/'
results_44249 <- 'results/psoriasis/magma/magma_output_50kb_window_snp-array_dsl2_only_human_ppi_44249/stable_consensus/data_splits/'
results <- 'results/psoriasis/magma/magma_output_50kb_window_snp-array_dsl2_only_human_ppi_44249/'

ppi <- read_tsv('data/biogrid_mv_physical_4.4.249_tab2_onlyhuman_9606_interactorsAB.tsv', col_types = 'cc') %>%
  select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
  graph_from_data_frame(directed = FALSE) %>%
  as_tbl_graph

methods <- c('Heinz','HotNet2', 'SigMod', 'cGWAS', 'MAGMABH', 'dmGWAS', 'LEAN', 'MAGMABHcGWAS')

method_palette <- c('Heinz' = '#73471f',
                    'HotNet2' = '#377eb8', 'SigMod' = '#f4888a', 'All' = 'gray50',
                    'Consensus' = 'black', 'MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'dmGWAS' = '#80cdc1', 'LEAN' = '#4d9221')

genesrc_palette <- c('MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'MAGMABH and cGWAS' = '#d3c783', 'Methods only' = '#ffffff')

# TOPOLOGY PANEL

plot_network <- function(net, selected_nodes, node_colors = NULL, label_genes = NULL) {
  
  graph <- activate(net, nodes) %>%
    filter(name %in% selected_nodes)
  
  graph <- graph %>%
    mutate(color = ifelse(name %in% names(node_colors), 
                          '#233142',
                          #node_colors[name], 
                          '#233142'))
  
  # ensure the graph is of class igraph and tbl_graph
  class(graph) <- c('igraph','tbl_graph')
  # graph <- as.igraph(graph)
  
  # get the ggnetwork data
  network_data <- ggnetwork(graph)
  
  # remove duplicate rows for labeling
  label_data <- network_data %>% filter(name %in% label_genes) %>% distinct(name, x, y)
  
  # plotting
  p <- ggplot(network_data, aes(x = x, y = y))
  
  # add edges only if xend and yend are present
  if("xend" %in% colnames(network_data) && "yend" %in% colnames(network_data)) {
    p <- p + geom_edges(aes(xend = xend, yend = yend), color = '#3282b8')
  }
  
  # add nodes, labels, and other styling
  p <- p +
    geom_nodes(aes(color = I(color)), size = 5) +
    geom_label_repel(
      data = label_data, 
      aes(label = name), 
      size = 3,
      box.padding = 0.2,
      point.padding = 0.3,
      segment.color = 'grey50'
    ) +
    ggtitle(' ') +
    theme_blank() + 
    theme_transparent + 
    theme(title = element_text(hjust = 1, face = 'bold'), legend.position="none")
  
  return(p)
}

# compilation of genes selected by all the methods across the different datasets.
heinz_genes <- c()
hotnet2_genes <- c()
sigmod_genes <- c()
dmgwas_genes <- c()
lean_genes <- c()

nsplits <- 5 # Number of splits applied to the data.
nmethods <- 5 # Number of biological network methods used

for (j in 1:nsplits) {
  # Read genes selected by Heinz
  file_name <- paste0(results_44249, "data_split_", j, "/heinz/selected_genes_split_", j, ".heinz.txt")
  heinz_genes <- unique(c(heinz_genes, read_tsv(file_name, col_types = 'c')$gene))
  
  # Read genes selected by HotNet2
  file_name <- paste0(results_44249, "data_split_", j, "/hotnet2/selected_genes_split_", j, ".hotnet2.tsv")
  hotnet2_genes <- unique(c(hotnet2_genes, read_tsv(file_name, col_types = 'cc')$gene))
  
  # Read genes selected by SigMod
  file_name <- paste0(results_44249, "data_split_", j, "/sigmod/selected_genes_split_", j, ".sigmod.txt")
  sigmod_genes <- unique(c(sigmod_genes, read_tsv(file_name, col_types = 'c')$gene))
  
  # Read genes selected by DMGWAS
  file_name <- paste0(results_44249, "data_split_", j, "/dmgwas/selected_genes_split_", j, ".dmgwas.txt")
  dmgwas_genes <- unique(c(dmgwas_genes, read_tsv(file_name, col_types = 'c')$gene))
  
  # Read genes selected by Lean
  file_name <- paste0(results_44249, "data_split_", j, "/lean/selected_genes_split_", j, ".lean.txt")
  lean_genes <- unique(c(lean_genes, read_tsv(file_name, col_types = 'c')$gene))
}


# powersets
gene_lists_subsets <- c()
gene_lists_src <- list()

for (j in 1:nsplits) {
  file_name <- paste0(results_44249,"data_split_", j, "/heinz/selected_genes_split_", j, ".heinz.txt")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'c')$gene

  file_name <- paste0(results_44249,"data_split_", j, "/hotnet2/selected_genes_split_", j, ".hotnet2.tsv")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'cc')$gene

  file_name <- paste0(results_44249,"data_split_", j, "/sigmod/selected_genes_split_", j, ".sigmod.txt")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'c')$gene
  
  file_name <- paste0(results_44249,"data_split_", j, "/dmgwas/selected_genes_split_", j, ".dmgwas.txt")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'c')$gene

  file_name <- paste0(results_44249,"data_split_", j, "/lean/selected_genes_split_", j, ".lean.txt")
  gene_lists_src[[length(gene_lists_src) + 1]] <- read_tsv(file_name, col_types = 'c')$gene
}

for (gene_list in gene_lists_src) {
  for (gene in gene_list) {
    if (gene %in% names(gene_lists_subsets)) {
      gene_lists_subsets[gene] <- gene_lists_subsets[gene] + 1
    } else {
      gene_lists_subsets[gene] <- 1
    }
  }
}

# Classical GWAS
cgwas_genes <- read_tsv(paste0(results_gwas,'pso.genes.annot_gene_ids_names.tsv'), col_types = 'c')$Gene %>% unique()
cgwas_plt <- plot_network(ppi, cgwas_genes)
cgwas_plt_square <- plot_network(ppi, cgwas_genes) + coord_fixed(ratio = 1)

# Magma output after Benjamini Hochberg correction
magmabh <- read_tsv(paste0(results,'pso.scores.genes.out_converted_pruned_bhcorrection.tsv'), 
                    col_types = 'iiciiiiiidnn')$Gene %>% unique()
magmabh_plt <- plot_network(ppi, magmabh)
magmabh_plt_square <- plot_network(ppi, magmabh) + coord_fixed(ratio = 1)

# Magma output intersected with classical gwas annotation output
magmabh_cgwas <- intersect(magmabh, cgwas_genes)
magmabh_cgwas_plt <- plot_network(ppi, magmabh_cgwas)
magmabh_cgwas_plt_square <- plot_network(ppi, magmabh_cgwas) + coord_fixed(ratio = 1)

# Calculate the number of genes appearing in at least 1, 2, 3, ..., nmethods * nsplits lists
max_occurrences <- nmethods * nsplits

# Loop through the process n times
for (i in seq(1, max_occurrences)) {
  # Genes in at least i solutions
  gene_set <- names(gene_lists_subsets[gene_lists_subsets >= i])
}

# Genes in at least 7 solutions

gene_set7 <- names(gene_lists_subsets[gene_lists_subsets >= 7])
gene_set7_plt <- plot_network(ppi, gene_set7)

# Gene set 7

gset7 <- activate(ppi, nodes) %>%
  filter(name %in% gene_set7) %>%
  mutate(gset7 = TRUE)
class(gset7) <- c('igraph','tbl_graph')

# cGWAS

g_cgwas <- activate(ppi, nodes) %>%
  filter(name %in% cgwas_genes) %>%
  mutate(cgwas_genes = TRUE)
class(g_cgwas) <- c('igraph','tbl_graph')

# MAGMA_BH

g_magmabh <- activate(ppi, nodes) %>%
  filter(name %in% magmabh) %>%
  mutate(magmabh = TRUE)
class(g_magmabh) <- c('igraph','tbl_graph')

# MAGMA_BH and cGWAS

g_magmabh_cgwas <- activate(ppi, nodes) %>%
  filter(name %in% magmabh_cgwas) %>%
  mutate(magmabh_cgwas = TRUE)
class(g_magmabh_cgwas) <- c('igraph','tbl_graph')

# Consensus

options(repr.plot.width=25, repr.plot.height=15)

consensus_gg <- function(gset, g_cgwas, g_magmabh, g_magmabh_cgwas, 
                         cgwas_genes, magmabh, magmabh_cgwas) {
  
  # gset <- gset %>% to_undirected
  
  gset <- gset %>% graph_join(g_magmabh_cgwas, by = c('name')) %>% 
    graph_join(g_cgwas, by = c('name')) %>%
    graph_join(g_magmabh, by = c('name')) %>%
    to_undirected %>%
    mutate(magmabh_cgwas = ifelse(is.na(magmabh_cgwas), FALSE, magmabh_cgwas),
           cgwas_genes = ifelse(is.na(cgwas_genes), FALSE, cgwas_genes),
           magmabh = ifelse(is.na(magmabh), FALSE, magmabh)) %>%
    filter(name %in% V(gset)$name)

  # output nodes selected by at least by 7 solutions
  as_tibble(gset) %>%
  rename(gene = name) %>%
  write_tsv(paste0(results_44249,'../stable_consensus_gset7.tsv'))
  
  set.seed(42)   # fixes the network layout for reproductibility
  gset %>%
  ggnetwork %>%
  mutate(name = as.character(name))
}

consensus <- consensus_gg(gset7, g_cgwas, g_magmabh, g_magmabh_cgwas, 
                          cgwas_genes, magmabh, magmabh_cgwas)

# output edge coordinates for each node in the stable consensus
consensus %>% select(x, y, name, xend, yend) %>%
  write_tsv(paste0(results_44249,'../connections_stable_consensus_gset7.tsv'))

nodes <- mutate(consensus,
                Heinz = as.numeric(name %in% heinz_genes),
                HotNet2 = as.numeric(name %in% hotnet2_genes),
                SigMod = as.numeric(name %in% sigmod_genes),
                dmGWAS = as.numeric(name %in% dmgwas_genes),
                LEAN = as.numeric(name %in% lean_genes),
                cGWAS = as.numeric(cgwas_genes),
                MAGMABH = as.numeric(magmabh),
                MAGMABHcGWAS = as.numeric(magmabh_cgwas)
                ) %>%
  filter(xend == x & yend == y) %>% 
  select(x, y, name, all_of(methods)) %>%
  unique

nodes <- nodes %>%
  mutate(label_group = case_when(
    MAGMABHcGWAS == 1 ~ 'MAGMABH and cGWAS',
    MAGMABH == 1 ~ 'MAGMABH',
    cGWAS == 1 ~ 'cGWAS',
    TRUE ~ 'Methods only')
  )

edges <- filter(consensus, xend != x | yend != y)

# CONSENSUS PIE

consensus_pie_methods <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.007), cols = setdiff(methods, c('MAGMABH', 'cGWAS', 'MAGMABHcGWAS'))) +
  scale_fill_manual(name = "Methods", values = method_palette) +
  ggnewscale::new_scale_fill() +
  geom_label_repel(data = nodes, 
                   aes(x = x, y = y, label = name, fill = label_group), 
                   size = 3, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 100) +
  scale_fill_manual(name = "Gene source", values = genesrc_palette) +
  coord_fixed() +
  theme_blank() +
  guides(color = "none", fill = guide_legend(
    override.aes = aes(label = "  ")
  )) +
  theme(legend.position = 'left',
        legend.text = element_text(size = 20, vjust = .2),
        legend.title = element_text(size = 22, vjust = 1)) +
  theme_transparent

consensus_pie_methods


ggsave(paste0(results_44249,'../figures/stable_consensus_methods_gset7.pdf'), consensus_pie_methods, width=14, height=10, limitsize = FALSE, bg = "transparent")
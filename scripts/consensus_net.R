# =============================================================================
# consensus_net.R  -  Consensus network of the genes selected by network methods
# =============================================================================
# Combines the modules returned by five network methods (Heinz, HotNet2, SigMod,
# dmGWAS and LEAN; Step 4 of gwas-bionets, bionets.nf) into a consensus network
# and draws it. Similarly, one can add more methods as needed to compare it in
# the consensus network. Two more gene sources are shown for comparison and annotation:
# genes annotated to the classical GWAS hits (cGWAS) and genes significant in
# MAGMA after Benjamini-Hochberg correction (MAGMABH).
# Author: Giann Karlo Aguirre-Samboní
#
# The general procedure of this script is as follows:
#   1. Loads a reference protein-protein interaction (PPI) network (e.g. BioGRID).
#   2. Reads the genes selected by each method (run labelled 'split_1', stored
#      directly in the results folder) plus the cGWAS and MAGMABH gene lists.
#   3. Builds the consensus: the genes selected by at least 3 of the 5 methods
#      (num_methods >= 3) together with the PPI edges among them. cGWAS and
#      MAGMABH are annotations and do not count towards the 3 methods.
#   4. Writes the consensus gene table and the node layout, and saves a figure in
#      which each node is a pie with one slice per method that selected the gene;
#      the label fill shows whether the gene is also in MAGMABH and/or cGWAS.
#
# Inputs (edit the paths below; '<prefix>' is the basename used in the pipeline)
#   PPI file       Tab-separated file whose first two columns are the official gene
#                  symbols of each interaction ('Official Symbol Interactor A' and
#                  'Official Symbol Interactor B', as in BioGRID tab2 files).
#   results_gwas   Folder with '<prefix>.genes.annot_gene_ids_names.tsv': table with a
#                  'Gene' column listing the genes annotated to the classical GWAS
#                  hits.
#   results_44249  Step 4 results folder (--out of bionets.nf; the name refers to the
#                  PPI release used). It must contain, for each method,
#                    heinz/selected_genes_split_1.heinz.txt
#                    hotnet2/selected_genes_split_1.hotnet2.tsv  (two columns)
#                    sigmod/selected_genes_split_1.sigmod.txt
#                    dmgwas/selected_genes_split_1.dmgwas.txt
#                    lean/selected_genes_split_1.lean.txt
#                  (each with a 'gene' column) and
#                    <prefix>.scores.genes.out_converted_pruned_bhcorrection.tsv
#                  MAGMA gene results after BH correction and pruning (deduplication): 12 columns, the
#                  third one being 'Gene' (character).
#
# Outputs (written to results_44249)
#   consensus.tsv              Genes of the consensus with one TRUE/FALSE column per
#                              source (heinz, hotnet2, sigmod, dmgwas, lean, magmabh_and_cgwas,
#                              cgwas_genes, magmabh) and num_methods.
#   connections_consensus.tsv  Layout coordinates of the nodes and edges (x, y, name, xend, yend).
#   figures/consensus_pie_magmabh_methods_cgwas.pdf  Figure (the 'figures' folder must exist, o/w 
#   modify the path).
#
# Notes
#   - The per-source plots (heinz_plt, hotnet2_plt, ...) are built for individual inspection if needed.
 
# Requirements: R packages igraph, tidygraph, magrittr, ggrepel, tidyverse, cowplot, network,
#   sna, ggnetwork, scatterpie and ggnewscale.
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
library(ggnewscale)

# Transparent ggplot theme (no panel, plot or legend background, no grid lines), so the figures can be placed on any background.
theme_transparent <- theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
                           plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                           panel.grid.major = element_blank(), # get rid of major grid
                           panel.grid.minor = element_blank(), # get rid of minor grid
                           legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
                           legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
)

theme_set(theme_cowplot() + theme_transparent)


results_gwas <- 'results/psoriasis/magma/magma_output_50kb_window_snp-array_significant_snps/'
results_44249 <- 'results/psoriasis/magma/magma_output_50kb_window_snp-array_dsl2_only_human_ppi_44249/'

ppi <- read_tsv('data/biogrid_mv_physical_4.4.249_tab2_onlyhuman_9606_interactorsAB.tsv', col_types = 'cc') %>%
  select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
  graph_from_data_frame(directed = FALSE) %>%
  as_tbl_graph

methods <- c('Heinz','HotNet2', 'SigMod', 'cGWAS', 'MAGMABH', 'dmGWAS', 'LEAN', 'MAGMABHcGWAS')

method_palette <- c('Heinz' = '#73471f',
                    'HotNet2' = '#377eb8', 'SigMod' = '#f4888a', 'All' = 'gray50',
                    'Consensus' = 'black', 'MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'dmGWAS' = '#80cdc1', 'LEAN' = '#4d9221')

genesrc_palette <- c('MAGMABH' = '#abdda4', 'cGWAS' = '#fdae61', 'MAGMABH and cGWAS' = '#d3c783', 'Methods only' = '#ffffff')

plot_network <- function(net, selected_nodes) {
  
  graph <- activate(net, nodes) %>%
    filter(name %in% selected_nodes)
  class(graph) <- c('igraph','tbl_graph')
  
  ggnetwork(graph) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = '#3282b8') +
    geom_nodes(color = '#233142', size = 4) +
    #geom_label(aes(label = name),
    #                 size = 2) +
    ggtitle(' ') +
    theme_blank() + theme_transparent + 
    theme(title = element_text(hjust = 1, face = 'bold'))
}

# heinz
heinz <- read_tsv(paste0(results_44249,'heinz/selected_genes_split_1.heinz.txt'), col_types = 'c')$gene
heinz_plt <- plot_network(ppi, heinz) 
heinz_plt_square <- heinz_plt + coord_fixed(ratio = 1) # Add coord_fixed(ratio = 1) to make the plot square

# hotnet2
hotnet2 <- read_tsv(paste0(results_44249,'hotnet2/selected_genes_split_1.hotnet2.tsv'), col_types = 'cc')$gene
hotnet2_plt <- plot_network(ppi, hotnet2)
hotnet2_plt_square <- hotnet2_plt + coord_fixed(ratio = 1)

# sigmod
sigmod <- read_tsv(paste0(results_44249,'sigmod/selected_genes_split_1.sigmod.txt'), col_types = 'c')$gene
sigmod_plt <- plot_network(ppi, sigmod)
sigmod_plt_square <- sigmod_plt + coord_fixed(ratio = 1)

# dmgwas
dmgwas <- read_tsv(paste0(results_44249,'dmgwas/selected_genes_split_1.dmgwas.txt'), col_types = 'c')$gene
dmgwas_plt <- plot_network(ppi, dmgwas)
dmgwas_plt_square <- dmgwas_plt + coord_fixed(ratio = 1)

# lean
lean <- read_tsv(paste0(results_44249,'lean/selected_genes_split_1.lean.txt'), col_types = 'c')$gene
lean_plt <- plot_network(ppi, lean)
lean_plt_square <- lean_plt + coord_fixed(ratio = 1)

# Classical GWAS
cgwas_genes <- read_tsv(paste0(results_gwas,'pso.genes.annot_gene_ids_names.tsv'), col_types = 'c')$Gene %>% unique()
cgwas_plt <- plot_network(ppi, cgwas_genes)
cgwas_plt_square <- plot_network(ppi, cgwas_genes) + coord_fixed(ratio = 1)

# Magma output after Benjamini Hochberg correction
magmabh <- read_tsv(paste0(results_44249,'pso.scores.genes.out_converted_pruned_bhcorrection.tsv'), 
                    col_types = 'iiciiiiiidnn')$Gene %>% unique()
magmabh_plt <- plot_network(ppi, magmabh)
magmabh_plt_square <- plot_network(ppi, magmabh) + coord_fixed(ratio = 1)

# Magma output together with classical gwas annotation output
# magmabh_cgwas <- union(magmabh, cgwas_genes)
# magmabh_cgwas_plt <- plot_network(ppi, magmabh_cgwas)
# magmabh_cgwas_plt_square <- plot_network(ppi, magmabh_cgwas) + coord_fixed(ratio = 1)

# Magma output intersected with classical gwas annotation output
magmabh_and_cgwas <- intersect(magmabh, cgwas_genes)
magmabh_and_cgwas_plt <- plot_network(ppi, magmabh_and_cgwas)
magmabh_and_cgwas_plt_square <- plot_network(ppi, magmabh_and_cgwas) + coord_fixed(ratio = 1)

# Heinz

g_heinz <- activate(ppi, nodes) %>%
  filter(name %in% heinz) %>%
  mutate(heinz = TRUE)
class(g_heinz) <- c('igraph','tbl_graph')

# HotNet2

g_hotnet2 <- activate(ppi, nodes) %>%
  filter(name %in% hotnet2) %>%
  mutate(hotnet2 = TRUE)
class(g_hotnet2) <- c('igraph','tbl_graph')

# SigMod

g_sigmod <- activate(ppi, nodes) %>%
  filter(name %in% sigmod) %>%
  mutate(sigmod = TRUE)
class(g_sigmod) <- c('igraph','tbl_graph')

# dmGWAS

g_dmgwas <- activate(ppi, nodes) %>%
  filter(name %in% dmgwas) %>%
  mutate(dmgwas = TRUE)
class(g_dmgwas) <- c('igraph','tbl_graph')

# LEAN

g_lean <- activate(ppi, nodes) %>%
  filter(name %in% lean) %>%
  mutate(lean = TRUE)
class(g_lean) <- c('igraph','tbl_graph')

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

# MAGMA_BH + cGWAS

# g_magmabh_cgwas <- activate(ppi, nodes) %>%
#   filter(name %in% magmabh_cgwas) %>%
#   mutate(magmabh_cgwas = TRUE)
# class(g_magmabh_cgwas) <- c('igraph','tbl_graph')

# MAGMA_BH and cGWAS

g_magmabh_and_cgwas <- activate(ppi, nodes) %>%
  filter(name %in% magmabh_and_cgwas) %>%
  mutate(magmabh_and_cgwas = TRUE)
class(g_magmabh_and_cgwas) <- c('igraph','tbl_graph')

# Consensus
# join all sources, set missing flags to FALSE, count how many of the five network methods
# selected each gene (num_methods) and keep genes selected by at least 3 of them.

options(repr.plot.width=25, repr.plot.height=15)
consensus <- graph_join(g_heinz, g_hotnet2, by = c('name')) %>%
  graph_join(g_sigmod, by = c('name')) %>%
  graph_join(g_dmgwas, by = c('name')) %>%
  graph_join(g_lean, by = c('name')) %>%
  graph_join(g_magmabh_and_cgwas, by = c('name')) %>%
  graph_join(g_cgwas, by = c('name')) %>%
  graph_join(g_magmabh, by = c('name')) %>%
  to_undirected %>%
  mutate(heinz = ifelse(is.na(heinz), FALSE, heinz),
         hotnet2 = ifelse(is.na(hotnet2), FALSE, hotnet2),
         sigmod = ifelse(is.na(sigmod), FALSE, sigmod),
         dmgwas = ifelse(is.na(dmgwas), FALSE, dmgwas),
         lean = ifelse(is.na(lean), FALSE, lean),
         magmabh_and_cgwas = ifelse(is.na(magmabh_and_cgwas), FALSE, magmabh_and_cgwas),
         cgwas_genes = ifelse(is.na(cgwas_genes), FALSE, cgwas_genes),
         magmabh = ifelse(is.na(magmabh), FALSE, magmabh),
         num_methods = rowSums(cbind(heinz, hotnet2, sigmod, dmgwas, lean))) %>%
  filter(num_methods >= 3)

# output nodes selected by different methods in TSV format
as_tibble(consensus) %>%
  rename(gene = name) %>%
  write_tsv(paste0(results_44249,'consensus.tsv'))

set.seed(42)   # fixes the network layout for reproductibility
consensus_gg <- consensus %>%
  ggnetwork %>%
  mutate(name = as.character(name))

# output edge coordinates for each node in the consensus
consensus_gg %>% select(x, y, name, xend, yend) %>%
  write_tsv(paste0(results_44249,'connections_consensus.tsv'))

# Label group (fill of the gene labels): MAGMABH and cGWAS > MAGMABH > cGWAS > methods only.
nodes <- mutate(consensus_gg, Heinz = as.numeric(heinz),
                HotNet2 = as.numeric(hotnet2),
                SigMod = as.numeric(sigmod),
                dmGWAS = as.numeric(dmgwas),
                LEAN = as.numeric(lean),
                cGWAS = as.numeric(cgwas_genes),
                MAGMABH = as.numeric(magmabh),
                MAGMABHcGWAS = as.numeric(magmabh_and_cgwas)
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

edges <- filter(consensus_gg, xend != x | yend != y)

consensus_pie_methods <- ggplot() +
  geom_edges(data = edges, aes(x = x, y = y, xend = xend, yend = yend), 
             color = '#393e46', linewidth = 1) +
  geom_scatterpie(data = nodes, aes(x = x, y = y, r=0.007), cols = setdiff(methods, c('MAGMABH', 'cGWAS', 'MAGMABHcGWAS'))) +
  scale_fill_manual(name = "Methods", values = method_palette) +
  ggnewscale::new_scale_fill() +
  geom_label_repel(data = nodes, 
                   aes(x = x, y = y, label = name, fill = label_group),
                   size = 3, point.padding = 1, max.time = 5, max.iter = 20000, 
                   max.overlaps = 50) +
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

ggsave(paste0(results_44249,'figures/consensus_pie_magmabh_methods_cgwas.pdf'), consensus_pie_methods, width=14, height=10, limitsize = FALSE, bg = "transparent")
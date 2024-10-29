#!/bin/sh
# Splitting scores
#
#
# Script: Data Preprocessing and Analysis Pipeline.
# Description: This script performs the following steps:
#   1. Splitting scores (or gene p-values)
#   2. Biological network construction
# Usage: ./bionets_construction_from_scores.sh
# Author: Giann Karlo Aguirre-Samboní
# Date: 20/10/2024


# Step 1: Splitting scores
# This step splits the list of gene scores (or gene p-values) into `k` parts that will be the input
# each one for the different algorithms.
# scrsample: path to where the scores are located.
# k: number of times scores will be split, it must be an integer gretar than 1.
# base_out_dir: address to save files, at least an output filename must be given.
# profile: nextflow variable to denote which setting use.
# dsl1: flag to restrict nextflow on domain-specific language 1 (DSL1), pipeline remains to be updated to DSL2.


scrsample=""
k=5
base_out_dir=""
profile="local"

bin/split_scores.nf \
  --scores "$scrsample" \
  --k "$k" \
  --out "$base_out_dir" \
  -dsl1 \
  -profile "$profile"

# Step 4: Construction of biological networks using different methods: HotNet2, SigMod and Heinz.

# net_ref: network of reference the methods will use for constructing subnetworks.
# net_results: path to folder to save results out of the methods.
# magma_scores: prefix of the filename where gene scores are stored and it is assumed to be located at net_results folder.
# This normally should be 'pso.scores.genes.out_converted' since magma_calc.nf output a file with this name.
# fdr: false discovery rate for use in Heinz method.
# lfdr_cutoff: local false discovery rate parameter of use for sparse scores of the HotNet2 method (cutoff if P-value >= 0.125).
# data_samp: boolean (1 or 0) to indicate whether the sampling (data splits) comes from the data or the gene scores (or p-values).
# sigmod: path to sigmod files with method's internal code.
# hotnet2: path to hotnet2 files with method's internal code.

net_ref=""
net_results=""
magma_scores="pso.scores.genes.out_converted"
fdr=0.5
lfdr_cutoff=0.125
data_samp=0
sigmod_path=""
hotnet2_path=""

# Biological network construction
bin/bionets.nf \
  --network "$net_ref" \
  --k "$k" \
  --d_samp "$data_samp" \
  --scores "$magma_scores" \
  --sigmod_path $sigmod_path \
  --hotnet2_path $hotnet2_path \
  --fdr "$fdr" \
  --lfdr "$lfdr_cutoff" \
  --out "$net_results" \
  -dsl1 \
  -profile "$profile"

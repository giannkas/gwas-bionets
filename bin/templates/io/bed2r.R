#!/usr/bin/env Rscript
# GENBED, VARBIM, SAMFAM
# gwas.RData

library(snpStats)

gwas <- read.plink("${GENBED}", "${VARBIM}", "${SAMFAM}")

save(gwas, file = 'gwas.RData')

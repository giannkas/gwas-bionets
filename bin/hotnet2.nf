#!/usr/bin/env nextflow

 params.out = "."
 params.hotnet2_path = null

 tab2 = file(params.tab2)
 magma = file(params.scores)
 params.i = 0
 params.d_samp = 1

 split = params.i > 0 && params.d_samp != 0 ? "_split_${params.i}" : params.i > 0 && params.d_samp == 0 ? "_chunk_${params.i}" : ""


 // conditional hotnet2 input to handle Docker or Dockerless execution
 docker_hotnet2 = '/gwas-bionets/hotnet2'

 if (params.hotnet2_path != null) {
     HOTNET2 = file(params.hotnet2_path)
 } else {
     HOTNET2 = docker_hotnet2
 }

 network_permutations = 100
 heat_permutations = 1000
 beta = 0.4
 params.lfdr_cutoff = 0.05

 process make_network {

     input:
         file TAB2 from tab2

     output:
         file 'node_index.tsv' into node_index
         file 'edge_list.tsv' into edge_list

     script:
     template 'io/tab2_2hotnet.R'

 }

 process sparse_scores {

     publishDir "$params.out", overwrite: true, mode: "copy"

     input:
     file SCORES from magma
     val CUTOFF from params.lfdr_cutoff
     val SPLIT from split

     output:
     file "scored_genes${SPLIT}.sparse.txt" into sparse_scores
     file "lfdr_plot${SPLIT}.pdf"

     """
#!/usr/bin/env Rscript

library(tidyverse)
library(twilight)
library(cowplot)

theme_set(theme_cowplot())

scores <- read_tsv('${SCORES}')

lfdr <- twilight(scores\$Pvalue, B=1000)
lfdr <- tibble(Gene = scores\$Gene[as.numeric(rownames(lfdr\$result))],
            magma_p = scores\$Pvalue[as.numeric(rownames(lfdr\$result))],
            lfdr = lfdr\$result\$fdr)

ggplot(lfdr, aes(x = magma_p, y = 1 - lfdr)) +
 geom_line() +
 geom_vline(xintercept = ${CUTOFF}, color = 'red') +
 labs(x = 'P-value', y = '1 - lFDR')
ggsave('lfdr_plot${SPLIT}.pdf', width=7, height=6)

lfdr %>%
 mutate(Pvalue = ifelse(magma_p < ${CUTOFF}, magma_p, 1)) %>%
 write_tsv('scored_genes${SPLIT}.sparse.txt')
     """

 }

 process magma2hotnet {

     input:
         file MAGMA from sparse_scores

     output:
         file 'scores.ht' into scores

     script:
     template 'io/magma2hotnet.R'

 }


 process make_h5_network {

     input:
         if (params.hotnet2_path != null) {
             file HOTNET2
         } else {
             val HOTNET2
         }
         file NODE_IDX from node_index
         file EDGE_LIST from edge_list
         val BETA from beta

     output:
         file "ppin_ppr_${BETA}.h5" into h5
         file 'permuted' into permutations

     """
     python2 ${HOTNET2}/makeNetworkFiles.py \
 --edgelist_file ${EDGE_LIST} \
 --gene_index_file ${NODE_IDX} \
 --network_name ppin \
 --prefix ppin \
 --beta ${BETA} \
 --cores -1 \
 --num_permutations ${network_permutations} \
 --output_dir .
     """

 }

 process make_heat_data {

     input:
         if (params.hotnet2_path != null) {
             file HOTNET2
         } else {
             val HOTNET2
         }
         file SCORES from scores

     output:
         file 'heat.json' into heat

     """
     python2 ${HOTNET2}/makeHeatFile.py \
 scores \
 --heat_file ${SCORES} \
 --output_file heat.json \
 --name gwas
     """

 }

 process hotnet2 {

     input:
         if (params.hotnet2_path != null) {
             file HOTNET2
         } else {
             val HOTNET2
         }
         file HEAT from heat
         file NETWORK from h5
         file PERMS from permutations
         val BETA from beta

     output:
         file 'consensus/subnetworks.tsv' into subnetworks

     """
     python2 ${HOTNET2}/HotNet2.py \
 --network_files ${NETWORK} \
 --permuted_network_path ${PERMS}/ppin_ppr_${BETA}_##NUM##.h5 \
 --heat_files ${HEAT} \
 --network_permutations ${network_permutations} \
 --heat_permutations ${heat_permutations} \
 --num_cores -1 \
 --output_directory .
     """

 }

 process process_output {

     publishDir "$params.out", overwrite: true, mode: "copy"

     input:
         file SUBNETWORKS from subnetworks
         val SPLIT from split

     output:
         file "selected_genes${SPLIT}.hotnet2.tsv" into genes_hotnet2

     """
 #!/usr/bin/env Rscript

 library(tidyverse)

 read_tsv('${SUBNETWORKS}', col_types = 'cc', comment = '#', col_names = F) %>%
     select(X1) %>%
     mutate(cluster = 1:n()) %>%
     separate_rows(X1, sep = ' ') %>%
     rename(gene = X1) %>%
     write_tsv('selected_genes${SPLIT}.hotnet2.tsv')
     """

 }


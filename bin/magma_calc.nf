#!/usr/bin/env nextflow

params.out = '.'
params.k = 5
params.window = 50

/////////////////////////////////////
//  MAGMA ANNOTATION
/////////////////////////////////////

process annotation_step {

  publishDir "${params.gene_annot}/data_splits/data_split_${task.index}", mode: "copy"

  input:
    val window from params.window
    val snplocpval from params.snploc_pval
    val geneloc from params.gene_loc
    val I from 1..params.k
    val K from params.k

  output:
    path "pso_split_${I}.genes.annot" into genesannot
    path "pso_split_${I}.log" into annotlog

  script:
    """
      /cluster/CBIO/data1/gaguirresamboni/ukb_gwas-tools/scripts/additional_tools/magma/magma \
        --annotate window=${window} \
        --snp-loc ${snplocpval}/data_splits/data_split_${I}/snp_association_plink_split_${I}.tsv \
        --gene-loc ${geneloc} \
        --out pso_split_${I} \
    """
}

/////////////////////////////////////
//  MAGMA GENE ANALYSIS
/////////////////////////////////////

process gene_analysis_step {

  publishDir "${params.out}/data_splits/data_split_${task.index}", mode: "copy"

  input:
    val I from 1..params.k
    val data from params.bfile
    val prx from params.prefix
    val snplocpval from params.snploc_pval
    path annot from genesannot

  output:
    path "pso.scores_split_${I}.genes.out" into genelist
    path "pso.scores_split_${I}.genes.raw" into generaw
    path "pso.scores_split_${I}.log" into genelog

  script:
    """
      lines=\$(grep -c '^' ${data}/data_splits/data_split_${I}/${prx}_split_${I}.fam)
      /cluster/CBIO/data1/gaguirresamboni/ukb_gwas-tools/scripts/additional_tools/magma/magma \
          --bfile ${data}/data_splits/data_split_${I}/${prx}_split_${I} \
          --pval ${snplocpval}/data_splits/data_split_${I}/snp_association_plink_split_${I}.tsv N=\$lines \
          --gene-annot ${annot} \
          --out pso.scores_split_${I}
    """

}

/////////////////////////////////////
//  FORMATTING GENE IDS FOR COMPATIBILITY
/////////////////////////////////////

process snps_pvalue_reformat {

  publishDir "${params.out}/data_splits/data_split_${task.index}", mode: "copy"

  input:
    val I from 1..params.k
    path genesmag from genelist

  output:
    path "pso.scores.genes.out_converted_split_${I}.tsv" into new_snppvalues

  script:
    """
      #!/usr/bin/env Rscript

      library(dplyr)
      magma <- readr::read_table("${genesmag}")
      gprofiler2::gconvert(magma$GENE, numeric_ns="ENTREZGENE_ACC") %>% 
        mutate(GENE = as.numeric(input)) %>% 
        select(GENE, name) %>% 
        left_join(magma, by = "GENE") %>%
        rename(Gene = name, Chr = CHR, nSNPs = NSNPS, Start = START, Stop = STOP, Test = ZSTAT, Pvalue = P) %>%
        readr::write_tsv("pso.scores.genes.out_converted_split_${I}.tsv")
    """

}

#!/usr/bin/env nextflow

log.info ""
log.info "========================================================"
log.info "|          [gwas-bionets] - magma_calc.nf          |"
log.info "========================================================"
log.info ""
log.info "### Annotation (to map SNPs onto genes) and ###"
log.info "### gene analysis steps (to compute gene p-values). ###"
log.info "### See MAGMAv1.10 user manual as reference. ###"
log.info ""
log.info ""
log.info ""
log.info "--------------------------------------------------------"
log.info "This program comes with NO WARRANTY"
log.info "It is free software, see LICENSE for details about"
log.info "redistribution and contribution."
log.info "--------------------------------------------------------"
log.info ""

// Help info ##########
params.help = null
if (params.help) {
    log.info ""
    log.info "Usage : magma_calc.nf --magma <binmagma> --window <window_size> --k <knumber> --snploc_pval <snps_file> --gene_loc <genes_file> \\"
    log.info "                      --bpfolder <genetic_data> --prefix <my_prefix> --plink <pversion> --gene_annot <annot_file> \\"
    log.info "                      --out <filename>"
    log.info ""
    log.info "  --magma       magma binary location, if not given its default value will be used, which is 'magma' called from PATH."
    log.info "  --window      number of base-pairs in kilobases (outside the transcription region) to be used as a range where"
    log.info "                SNPs may influence a gene."
    log.info "  --k           number of times data was split, it must be in accordance to the previous k for splitting the data."
    log.info "  --snploc_pval SNP location file, it can be the file from PLINK snp-pvalues or .bim file (see MAGMA manual for more information)."
    log.info "  --gene_loc    gene location file, relative to a particular human genome reference (eg. NCBI37.3.gene.loc)."
    log.info "  --bpfolder    path to parent folder where genotypes has been split and where 'data_splits' folder resides so it can "
    log.info "                search iteratively the genetic data input (bed, bim and fam files) for each data split (1, 2, ..., k)."
    log.info "                If there was no splitting process, then providing the parent folder will suffice. MAGMA works on PLINK1 format."
    log.info "  --gene_annot  path to parent folder where the corresponding annotation files has been created depending on the splits."
    log.info "                This folder should also point to the place where 'data_splits' folder is, normally different from that of the"
    log.info "                genotypes because it is part of the analysis outputs."
    log.info "  --prefix      prefix or basename for the genetic data input (eg. my_prefix.bed, my_prefix.bim, my_prefix.fam)."
    log.info "  --plink       version of PLINK to use (1 or 2), note that pgen, pvar and psam is for PLINK2, whereas"
    log.info "                bed, bim and fam corresponds to PLINK. Since MAGMA only works on PLINK format so this value must be 1."
    log.info "  --out         path/to/filename where output files will be saved."
    log.info ""
    log.info ""
    log.info ""
    log.info "Example : "
    log.info "split_data.nf \\"
    log.info "  --window            50 \\"
    log.info "  --k                 5 \\"
    log.info "  --snploc_pval       path/to/snp_locations \\"
    log.info "  --gene_loc          path/to/gene_locations \\"
    log.info "  --bpfolder          path/to/genotypes_folder \\"
    log.info "  --gene_annot        path/to/annotations_folder \\"
    log.info "  --prefix            my_prefix \\"
    log.info "  --plink             1 \\"
    log.info "  --out               path/to/my_output \\"
    log.info "  -dsl1 \\"
    log.info "  -profile            my_cluster \\"
    log.info ""

    exit 0
}



params.out = '.'
params.k = 1
params.window = 50
params.plink = 1
params.magma = "magma"

/////////////////////////////////////
//  MAGMA ANNOTATION
/////////////////////////////////////

process annotation_step {

  publishDir { params.out + (params.k > 1 ? '/data_splits/data_split_' + task.index : '') }, mode: "copy"

  input:
    val window from params.window
    val snplocpval from params.snploc_pval
    val geneloc from params.gene_loc
    val I from 1..params.k
    val K from params.k
    val magma from params.magma

  output:
    path {"pso" + (K > 1 ? "_split_" + I : '') + ".genes.annot" } into genesannot
    path {"pso" + (K > 1 ? "_split_" + I : '') + ".log" } into annotlog

  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""
    def dir_splits = K > 1 ? "/data_splits/data_split_${I}" : ""

    """
      ${magma} \
        --annotate window=${window} \
        --snp-loc ${snplocpval}${dir_splits}/snp_association_plink${split_suffix}.tsv \
        --gene-loc ${geneloc} \
        --out pso${split_suffix} \
    """
}

/////////////////////////////////////
//  MAGMA GENE ANALYSIS
/////////////////////////////////////

process gene_analysis_step {

  publishDir { params.out + (params.k > 1 ? '/data_splits/data_split_' + task.index : '') }, mode: "copy"

  input:
    val I from 1..params.k
    val data from params.bpfolder
    val prx from params.prefix
    val snplocpval from params.snploc_pval
    path annot from genesannot
    val K from params.k
    val V from params.plink
    val magma from params.magma

  output:
    path {"pso.scores" + (K > 1 ? "_split_" + I : '') + ".genes.out" } into genelist
    path {"pso.scores" + (K > 1 ? "_split_" + I : '') + ".genes.raw" } into generaw
    path {"pso.scores" + (K > 1 ? "_split_" + I : '') + ".log" } into genelog


  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""
    def dir_splits = K > 1 ? "/data_splits/data_split_${I}" : ""

    """
    if [ $V -eq 1 ]; then
      lines=\$(grep -c '^' ${data}${dir_splits}/${prx}${split_suffix}.fam)
    else
      lines=\$(grep -c '^' ${data}${dir_splits}/${prx}${split_suffix}.psam)
    fi

      ${magma} \
        --bfile ${data}${dir_splits}/${prx}${split_suffix} \
        --pval ${snplocpval}${dir_splits}/snp_association_plink${split_suffix}.tsv N=\$lines \
        --gene-annot ${annot} \
        --out pso.scores${split_suffix}
    """

}

//////////////////////////////////////////////
//  FORMATTING GENE IDS FOR COMPATIBILITY
//////////////////////////////////////////////

process snps_pvalue_reformat {

  publishDir { params.out + (params.k > 1 ? '/data_splits/data_split_' + task.index : '') }, mode: "copy"


  input:
    val I from 1..params.k
    path genesmag from genelist
    val K from params.k

  output:
    path {"pso.scores.genes.out_converted" + (K > 1 ? "_split_" + I : '') + ".tsv" } into new_snppvalues

  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""
    def dir_splits = K > 1 ? "/data_splits/data_split_${I}" : ""

    """
      #!/usr/bin/env Rscript

      library(dplyr)
      magma <- readr::read_table("${genesmag}")
      gprofiler2::gconvert(magma\$GENE, numeric_ns="ENTREZGENE_ACC") %>% 
        mutate(GENE = as.numeric(input)) %>% 
        select(GENE, name) %>% 
        left_join(magma, by = "GENE") %>%
        rename(Gene = name, Chr = CHR, nSNPs = NSNPS, Start = START, Stop = STOP, Test = ZSTAT, Pvalue = P) %>%
        readr::write_tsv("pso.scores.genes.out_converted${split_suffix}.tsv")
    """

}

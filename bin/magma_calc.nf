#!/usr/bin/env nextflow

// =============================================================================
// magma_calc.nf  -  Step 3 of gwas-bionets: MAGMA annotation and gene analysis
// =============================================================================
// Uses MAGMA (v1.10) to map SNPs onto genes and to turn SNP P-values into
// gene-level association statistics, for each data split. See the MAGMA manual.
//
// Workflow
//   1. annotation_step     magma --annotate window=<W> ...  -> '<prefix>[_split_<I>].genes.annot'
//   2. gene_analysis_step  magma --bfile ... --pval <snp P-values> N=<samples> --gene-annot ...
//                          -> '<prefix>.scores[_split_<I>].genes.out' (+ '.genes.raw', '.log').
//                          The genotypes are the LD reference; N is the number of samples
//                          (lines of the .fam file) of the split.
//   3. gene_ids_reformat   R script: converts the Entrez gene IDs of MAGMA's output into gene
//                          symbols (g:Profiler, R package gprofiler2: needs internet access),
//                          keeps one row per gene (the smallest P-value), drops unmapped or
//                          incomplete rows and renames the columns to Gene, Chr, nSNPs, Start,
//                          Stop, Test (z-statistic) and Pvalue
//                          -> '<prefix>.scores.genes.out_converted[_split_<I>].tsv'
//
// Parameters (command line, as --name value)
//   --magma        Path to the MAGMA binary (default 'magma', taken from the PATH).
//   --window       Window, in kilobases, added upstream and downstream of each gene when
//                  annotating SNPs (default 50).
//   --k            Number of splits; it must match the --k used in split_data.nf (default 1).
//   --snploc_pval  Folder with the SNP association files written by snps_pvalue.nf
//                  ('<prefix>_snp_association_plink[_split_<I>].tsv', in
//                  <folder>/stable_consensus/data_splits/data_split_<I>/ when k > 1).
//   --gene_loc     Gene location file for the genome build of the data (e.g. NCBI37.3.gene.loc).
//   --bpfolder     Parent folder of the PLINK 1 genotype files (bed/bim/fam); with k > 1 they
//                  are read from <bpfolder>/data_splits/data_split_<I>/<prefix>_split_<I>.
//   --prefix       Basename of the genetic data input.
//   --out          Output folder (default '.'); with k > 1 results go to
//                  <out>/stable_consensus/data_splits/data_split_<I>/.
//   --gene_annot, --plink   Accepted for compatibility but NOT read by this script:
//                  the annotation files come from the annotation step of the same run, and
//                  MAGMA always works on PLINK 1 files.
//   --help         Print the usage and exit.
//
// Requirements: MAGMA; R with dplyr, tidyr, readr and gprofiler2.

/////////////////////////////////////
//  MAGMA ANNOTATION
/////////////////////////////////////

// annotation_step: one task per split I (1..K).
//   window      --window, in kb
//   snplocpval  --snploc_pval (folder with the SNP association files)
//   geneloc     --gene_loc
//   magma       --magma, K, I, prefix: MAGMA binary, number of splits, split index, --prefix
//   emits       genesannot ('.genes.annot') and annotlog ('.log')
process annotation_step {

  publishDir { params.out + ((params.k as Integer) > 1 ? '/stable_consensus/data_splits/data_split_' + I : '') }, mode: "copy"

  input:
    val window
    path snplocpval
    path geneloc
    val magma
    val K
    each I
    val prefix

  output:
    path {"${prefix}" + (K > 1 ? "_split_" + I : '') + ".genes.annot" }, emit: genesannot
    path {"${prefix}" + (K > 1 ? "_split_" + I : '') + ".log" }, emit: annotlog

  script:
    def split_suffix = K > 1 ? "_split_${I}" : ""
    def dir_splits = K > 1 ? "/stable_consensus/data_splits/data_split_${I}" : ""

    """
      $magma \
        --annotate window=${window} \
        --snp-loc ${snplocpval}${dir_splits}/${prefix}_snp_association_plink${split_suffix}.tsv \
        --gene-loc ${geneloc} \
        --out ${prefix}${split_suffix}
    """
}

/////////////////////////////////////
//  MAGMA GENE ANALYSIS
/////////////////////////////////////

// gene_analysis_step: one task per split.
//   data        --bpfolder (genotypes, LD reference)
//   genesannot  annotation file(s) from annotation_step
//   snplocpval  folder with the SNP P-values; N = number of lines of the split's .fam file
//   emits       genesmag ('.genes.out'), generaw ('.genes.raw') and genelog ('.log')
process gene_analysis_step {

  publishDir { params.out + ((params.k as Integer) > 1 ? '/stable_consensus/data_splits/data_split_' + task.index : '') }, mode: "copy"

  input:
    path data
    path genesannot
    val prefix
    path snplocpval
    val magma
    val K

  output:
    path {"${prefix}.scores" + (K > 1 ? "_split_" + task.index : '') + ".genes.out" }, emit: genesmag
    path {"${prefix}.scores" + (K > 1 ? "_split_" + task.index : '') + ".genes.raw" }, emit: generaw
    path {"${prefix}.scores" + (K > 1 ? "_split_" + task.index : '') + ".log" }, emit: genelog


  script:
    def split_suffix = K > 1 ? "_split_${task.index}" : ""
    def dir_splits = K > 1 ? "/data_splits/data_split_${task.index}" : ""
    def scdir_splits = K > 1 ? "/stable_consensus/data_splits/data_split_${task.index}" : ""

    // if [ $vplink -eq 1 ]; then
    //   lines=\$(grep -c '^' ${data}${dir_splits}/${prefix}${split_suffix}.fam)
    // else
    //   lines=\$(grep -c '^' ${data}${dir_splits}/${prefix}${split_suffix}.psam)
    // fi

    """
      lines=\$(grep -c '^' ${data}${dir_splits}/${prefix}${split_suffix}.fam)

      $magma \
        --bfile ${data}${dir_splits}/${prefix}${split_suffix} \
        --pval ${snplocpval}${scdir_splits}/${prefix}_snp_association_plink${split_suffix}.tsv N=\$lines \
        --gene-annot ${genesannot} \
        --out ${prefix}.scores${split_suffix}
    """

}


//////////////////////////////////////////////
//  FORMATTING GENE IDS FOR COMPATIBILITY
//////////////////////////////////////////////

// gene_ids_reformat: converts Entrez IDs to gene symbols and tidies the gene table
// (see the header). emits new_geneids: '<prefix>.scores.genes.out_converted[_split_<n>].tsv'
process gene_ids_reformat {

  publishDir { params.out + ((params.k as Integer) > 1 ? '/stable_consensus/data_splits/data_split_' + task.index : '') }, mode: "copy"


  input:
    path genesmag
    val K
    val prefix

  output:
    path {"${prefix}.scores.genes.out_converted" + (K > 1 ? "_split_" + task.index : '') + ".tsv" }, emit: new_geneids

  script:
    def split_suffix = K > 1 ? "_split_${task.index}" : ""
    def dir_splits = K > 1 ? "/data_splits/data_split_${task.index}" : ""

    """
      #!/usr/bin/env Rscript

      library(dplyr)
      library(tidyr)
      magma <- readr::read_table("${genesmag}")
      gprofiler2::gconvert(magma\$GENE, numeric_ns="ENTREZGENE_ACC", mthreshold=1, filter_na=TRUE) %>% 
        mutate(GENE = as.numeric(input)) %>% 
        select(GENE, name) %>% 
        left_join(magma, by = "GENE") %>%
        rename(Gene = name, Chr = CHR, nSNPs = NSNPS, Start = START, Stop = STOP, Test = ZSTAT, Pvalue = P) %>%
        filter(Gene != "None") %>%
        group_by(Gene) %>%
        filter(Pvalue == min(Pvalue)) %>%
        ungroup() %>%
        distinct(Gene, .keep_all = TRUE) %>%
        drop_na() %>%
        readr::write_tsv("${prefix}.scores.genes.out_converted${split_suffix}.tsv")
    """

}

params.out = '.'             // --out: output folder
params.k = 1                 // --k: number of splits (same as in split_data.nf)
params.window = 50            // --window: SNP-to-gene window in kb (up- and downstream)
params.plink = 1             // --plink: not read by this script (MAGMA uses PLINK 1 files)
params.magma = "magma"       // --magma: path to the MAGMA binary
params.snploc_pval = ""       // --snploc_pval: folder with the SNP association files (Step 2 output)
params.gene_loc = ""          // --gene_loc: gene location file
params.bpfolder = ""         // --bpfolder: parent folder of the PLINK 1 genotype files
params.prefix = ""           // --prefix: basename of the genetic data input

// Help info ##########
params.help = null            // --help: print the usage and exit

// Workflow
workflow {
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

  if (params.help) {
    log.info "Usage : magma_calc.nf --magma <binmagma> --window <window_size> --k <knumber> --snploc_pval <snps_folder> --gene_loc <genes_file> \\"
    log.info "                      --bpfolder <genetic_data> --prefix <my_prefix> --out <filename>"
    log.info ""
    log.info "  --magma       magma binary location, if not given its default value will be used, which is 'magma' called from PATH."
    log.info "  --window      window size, in kilobases, added upstream and downstream of each gene (outside the transcription"
    log.info "                region) to be used as a range where SNPs may influence a gene (default 50)."
    log.info "  --k           number of times data was split, it must be in accordance to the previous k for splitting the data."
    log.info "  --snploc_pval folder with the SNP association files produced by snps_pvalue.nf (the --out of that step),"
    log.info "                '<prefix>_snp_association_plink[_split_<I>].tsv' (columns SNP, CHR, BP, P). With k > 1 they are"
    log.info "                searched in <folder>/stable_consensus/data_splits/data_split_<I>/."
    log.info "  --gene_loc    gene location file, relative to a particular human genome reference (eg. NCBI37.3.gene.loc)."
    log.info "                Its genome build must match that of the genotype data."
    log.info "  --bpfolder    path to parent folder where genotypes has been split and where 'data_splits' folder resides so it can"
    log.info "                search iteratively the genetic data input (bed, bim and fam files) for each data split (1, 2, ..., k)."
    log.info "                If there was no splitting process, then providing the parent folder will suffice. MAGMA works on PLINK1 format."
    log.info "  --prefix      prefix or basename for the genetic data input (eg. my_prefix.bed, my_prefix.bim, my_prefix.fam)."
    log.info "  --out         path/to/folder where output files will be saved. With k > 1 results go to"
    log.info "                <out>/stable_consensus/data_splits/data_split_<I>/."
    log.info "  --gene_annot  accepted for compatibility but not read: annotation files are taken from the annotation step of the same run."
    log.info "  --plink       accepted for compatibility but not read: MAGMA only works on PLINK 1 format (bed, bim and fam)."
    log.info "  --help        show this message and exit."
    log.info ""
    log.info "Example :"
    log.info "magma_calc.nf \\"
    log.info "  --magma             path/to/bin/magma \\"
    log.info "  --window            50 \\"
    log.info "  --k                 5 \\"
    log.info "  --snploc_pval       path/to/snp_results_folder \\"
    log.info "  --gene_loc          path/to/gene_locations \\"
    log.info "  --bpfolder          path/to/genotypes_folder \\"
    log.info "  --prefix            my_prefix \\"
    log.info "  --out               path/to/my_output \\"
    log.info "  -profile            my_cluster \\"

    exit 0
  }

  // Define inputs
  def data = params.bpfolder
  def window = params.window
  def snplocpval = file(params.snploc_pval)
  def geneloc = file(params.gene_loc)
  def magma = file(params.magma)
  def K = params.k as Integer
  def prefix = params.prefix

  // println "The value of magma is: ${magma}"
  annotation_step(
    window, snplocpval, geneloc, magma, K, 1..K, prefix)

  gene_analysis_step(
    data, annotation_step.out.genesannot, prefix, snplocpval, magma, K)

  gene_ids_reformat(gene_analysis_step.out.genesmag, K, prefix)

}
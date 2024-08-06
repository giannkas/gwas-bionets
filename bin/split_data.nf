#!/usr/bin/env nextflow

params.out = '.'
params.k = 5

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.baseName}.bim")
fam = file("${params.bfile}.fam")

/////////////////////////////////////
//  SPLIT PREPARATION
/////////////////////////////////////
process make_splits {

  publishDir "${params.out}/data_splits/data_split_${task.index}", overwrite: true, mode: "copy"

  input:
    file FAM from fam
    val I from 1..params.k
    val K from params.k

  output:
    path "${bed.baseName}_split_${I}.fam" into splits

  script:
  template 'genotypes/split_fams.sh'
  

}


/////////////////////////////////////
//  BIM AND BED GENERATION FROM FAM
/////////////////////////////////////

process bim_bed_creation {

  publishDir "${params.out}/data_splits/data_split_${task.index}", mode: "copy"

  input:
    val data from params.bfile
    val I from 1..params.k
    val K from params.k
    path fams from splits

  output:
    path "${bed.baseName}_split_${I}.bim" into bims
    path "${bed.baseName}_split_${I}.bed" into beds

  script:
    """
      plink --bfile ${data} --keep ${fams} --make-bed --out ${bed.baseName}_split_${I}
    """
}


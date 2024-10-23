#!/usr/bin/env nextflow

params.out = '.'
params.k = 1
params.scores = ""

tsv = file("${params.scores}")
filename = tsv.baseName

/////////////////////////////////////
//  CHUNK PREPARATION
/////////////////////////////////////
process make_chunks {

  publishDir { params.out + (params.k > 1 ? '/scores_chunks/scores_chunk_' + task.index : '') }, overwrite: true, mode: "copy"

  input:
    val SCORES from tsv
    val I from 1..params.k
    val K from params.k
    val prx from filename

  output:
    path { filename + (K > 1 ? "_chunk_" + I : '') + ".tsv" } into chunks

  script:
    //def chunk_suffix = K > 1 ? "_chunk_${I}" : ""

  // template 'genotypes/chunk_scores.sh'

  
  """
    #!/usr/bin/env Rscript

      library(tidyverse)
      scores <- read_tsv('${SCORES}')

      #define number of data frames to split into
      n <- ${K}
      if (n > 1) {
        chunk_suffix <- paste0("_chunk_", ${I})
      } else {
        chunk_suffix <- ""
      }
      nrows <- nrow(scores)

      #split data frame into n equal-sized data frames
      #chunks <- split(scores, factor(sort(rank(row.names(scores))%%n)))

      #loop through each chunk and save as tsv file
      #for (i in seq_along(chunks)) {

      # Calculate the range of rows to exclude in this chunk
      exclude_start <- ceiling((${I} - 1) * nrows / n) + 1
      exclude_end <- ceiling(${I} * nrows / n)

      # Create the chunk by excluding the specified rows
      chunk <- scores[-(exclude_start:exclude_end), ]

      #create filename based on chunk number
      outfile <- paste0("${prx}", chunk_suffix, ".tsv")
      
      #write.table(chunks[[i]], file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
      #write the chunk to a tsv file with row.names=FALSE to avoid saving row names
      write.table(chunk, file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      #}
  """
  

}


create_illumina_quality_frame <- function(
    input_path_transcriptome,
    first_lengths,
    second_lengths,
    first_quality,
    second_quality,
    alignment,
    sample_number,
    seed,
    output_path) {
  set.seed(seed)

  transcriptome <- data.table::fread(
    input_path_transcriptome
  )

  selected_names <- sample(transcriptome$qname, sample_number)

  length_first <- data.table::fread(
    first_lengths
  )

  length_second <- data.table::fread(
    second_lengths
  )

  quality_first <- data.table::fread(
    first_quality
  )

  quality_second <- data.table::fread(
    second_quality
  )
  
  alignment <- data.table::fread(
    alignment
  )


  transcriptome <- transcriptome[match(selected_names, transcriptome$qname), ]
  length_first <- length_first[match(selected_names, length_first$V1), ]
  length_second <- length_second[match(selected_names, length_second$V1), ]

  quality_first <- quality_first[match(selected_names, quality_first$V1), ]
  quality_second <- quality_second[match(selected_names, quality_second$V1), ]
  
  alignment <- alignment[match(selected_names, alignment$qname), ]

  data.frame(
    qname = selected_names,
    tlen = abs(transcriptome$tlen),
    length = length_first$V2 + length_second$V2,
    quality = apply(cbind(quality_first$V2, quality_second$V2), 1, mean),
    n_junctions = alignment$n_junctions,
    edit_distance = alignment$edit_distance
  ) %>% write_tsv(output_path)
}



log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
})

status <- create_illumina_quality_frame(
  input_path_transcriptome = snakemake@input[["input_path_transcriptome"]],
  first_lengths = snakemake@input[["first_lengths"]],
  second_lengths = snakemake@input[["second_lengths"]],
  first_quality = snakemake@input[["first_quality"]],
  second_quality = snakemake@input[["second_quality"]],
  alignment = snakemake@input[["alignment"]],
  sample_number = snakemake@params[["sample_number"]],
  seed = snakemake@params[["seed"]],
  output_path = snakemake@output[[1]]
)

sink()
sink()

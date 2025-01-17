log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

prepare_sqanti_protein <- function(query, reference, query_output_path, reference_output_path) {
  query <- import(query)
  query <- query[query$type == "CDS"]
  query$type <- "exon"
  query$phase <- NA

  export(query, query_output_path)

  reference <- import(reference)

  reference <- reference[reference$type == "CDS"]
  reference$type <- "exon"
  reference$phase <- NA

  export(reference, reference_output_path)
  return(0)
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
})

status <- prepare_sqanti_protein(
  query = snakemake@input[["query"]],
  reference = snakemake@input[["reference"]],
  query_output_path = snakemake@output[["query"]],
  reference_output_path = snakemake@output[["reference"]]
)

sessionInfo()

sink()
sink()

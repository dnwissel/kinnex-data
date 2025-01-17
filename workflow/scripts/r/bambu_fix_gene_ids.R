bambu_fix_gene_ids <- function(input_path, output_path) {
  transcriptome <- import(input_path)
  gene_length <- sapply(strsplit(transcriptome$gene_id, ";"), length)
  
  transcriptome$gene_id[gene_length != 1] <- sapply(strsplit(transcriptome$gene_id[gene_length != 1], ";"), function(x) trimws(x[[2]]))
  transcriptome <- transcriptome[strand(transcriptome) != "*"]
  export(transcriptome, output_path)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages(library(rtracklayer))

status <- bambu_fix_gene_ids(
  input_path = snakemake@input[["input_path"]],
  output_path = snakemake@output[["output_path"]]
)

sink()
sink()

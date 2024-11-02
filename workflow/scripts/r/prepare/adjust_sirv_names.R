adjust_sirv_names <- function(input_path, genome_name, transcriptome_name, output_path_transcriptome, output_path_genome) {
  genome <- readDNAStringSet(paste0(input_path, "/", genome_name))
  transcriptome <- import(paste0(input_path, "/", transcriptome_name))

  transcriptome$transcript_id[which(transcriptome$transcript_id %in% names(genome))] <- paste0(transcriptome$transcript_id[which(transcriptome$transcript_id %in% names(genome))], "_transcript")
  transcriptome$gene_id[which(transcriptome$gene_id %in% names(genome))] <- paste0(transcriptome$gene_id[which(transcriptome$gene_id %in% names(genome))], "_gene")
  export(transcriptome, output_path_transcriptome)
  writeXStringSet(genome, output_path_genome)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
})

status <- adjust_sirv_names(
  input_path = snakemake@input[[1]],
  genome_name = snakemake@params[["genome_name"]],
  transcriptome_name = snakemake@params[["transcriptome_name"]],
  output_path_transcriptome = snakemake@output[["output_path_transcriptome"]],
  output_path_genome = snakemake@output[["output_path_genome"]]
)

sink()
sink()

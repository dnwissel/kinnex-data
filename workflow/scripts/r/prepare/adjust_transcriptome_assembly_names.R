adjust_transcriptome_assembly_names <- function(genome, transcriptome, output_path) {
  genome <- readDNAStringSet(genome)
  transcriptome <- import(transcriptome)
  new_levels <- levels(seqnames(transcriptome))
  for (contig in unique(seqnames(transcriptome))[!(unique(seqnames(transcriptome)) %in% c(paste0("chr", c(1:22, "X", "Y", "M"))))]) {
    new_levels <- c(new_levels, grep(contig, sapply(strsplit(names(genome), " "), function(x) x[[1]]), value = TRUE))
  }
  seqlevels(transcriptome) <- new_levels
  for (contig in unique(seqnames(transcriptome))[!(unique(seqnames(transcriptome)) %in% c(paste0("chr", c(1:22, "X", "Y", "M"))))]) {
    seqnames(transcriptome[seqnames(transcriptome) == contig]) <- factor(grep(contig, sapply(strsplit(names(genome), " "), function(x) x[[1]]), value = TRUE), levels = seqlevels(transcriptome))
  }

  seqlevels(transcriptome) <- levels(droplevels(seqnames(transcriptome)))
  export(transcriptome, output_path)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
})

status <- adjust_transcriptome_assembly_names(
  genome = snakemake@input[["genome"]],
  transcriptome = snakemake@input[["transcriptome"]],
  output_path = snakemake@output[[1]]
)

sink()
sink()

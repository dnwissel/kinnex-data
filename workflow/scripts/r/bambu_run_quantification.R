bambu_run_quantification <- function(reads, annotations, genome, discovery, output_directory, ncore, NDR, quantification, lowMemory = TRUE) {
  bambuAnnotations <- prepareAnnotations(annotations)
  bambu_result <- bambu::bambu(
    reads = reads, annotations = bambuAnnotations,
    genome = genome, discovery = discovery,
    ncore = ncore, lowMemory = lowMemory, NDR = NDR, quant = quantification
  )
  bambu::writeBambuOutput(bambu_result, path = output_directory)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages(library(bambu))

if (snakemake@params[["sirv"]] == "sirv") {
  annotation <- snakemake@input[["sirv_transcriptome"]]
  reference <- snakemake@input[["sirv_genome"]]
} else {
  annotation <- snakemake@input[["gencode_transcriptome"]]
  reference <- snakemake@input[["gencode_genome"]]
}

status <- bambu_run_quantification(
  reads = unlist(strsplit(snakemake@input[["reads"]], ",")),
  annotations = annotation,
  genome = reference,
  discovery = FALSE,
  output_directory = snakemake@params[["outdir"]],
  ncore = snakemake@threads,
  lowMemory = TRUE,
  NDR = NULL,
  quantification = TRUE
)

sink()
sink()

run_bambu <- function(reads, annotations, genome, output_path, ncore, NDR, quantification, lowMemory = TRUE) {
  bambu_result <- bambu::bambu(
    reads = reads, annotations = annotations,
    genome = genome, discovery = TRUE,
    ncore = ncore, lowMemory = lowMemory, NDR = NDR, quant = FALSE,
    opt.discovery = list(remove.subsetTx = FALSE)
  )
  writeToGTF(bambu_result, file = output_path)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(bambu))

print(snakemake@input[["reads"]])
sessionInfo()
status <- run_bambu(
  reads = snakemake@input[["reads"]],
  annotations = NULL,
  genome = snakemake@input[["genome"]],
  output_path = snakemake@output[[1]],
  ncore = snakemake@threads,
  lowMemory = TRUE,
  NDR = 1.0
)

sink()
sink()

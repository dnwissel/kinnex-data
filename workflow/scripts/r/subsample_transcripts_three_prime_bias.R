log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

subsample_transcripts_three_prime_bias <- function(
    transcriptome,
    subsampled_transcriptome_path,
    n_subsample,
    seed) {
  set.seed(as.integer(seed))
  transcriptome <- import(transcriptome)
  transcript_ids <- unique(transcriptome$transcript_id)
  sampled_transcripts <- sample(x = transcript_ids, replace = FALSE, size = n_subsample)
  subsampled_transcriptome <- transcriptome[transcriptome$transcript_id %in% sampled_transcripts]
  export(
    subsampled_transcriptome,
    subsampled_transcriptome_path
  )
  return(0)
}

suppressPackageStartupMessages(library(rtracklayer))

status <- subsample_transcripts_three_prime_bias(
  transcriptome = snakemake@input[["transcriptome"]],
  subsampled_transcriptome_path = snakemake@output[["subsampled_transcriptome"]],
  n_subsample = snakemake@params[["n_subsample"]],
  seed = snakemake@params[["seed"]]
)

sessionInfo()

sink()
sink()

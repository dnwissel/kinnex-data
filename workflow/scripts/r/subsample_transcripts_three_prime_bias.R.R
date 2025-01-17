log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

subsample_bam <- function(
    input_name,
    output_name,
    n_unique_output_reads,
    read_name_sam,
    seed,
    is_transcriptome) {
  set.seed(as.integer(seed))
  input_reads <- unname(unlist(scanBam(file = input_name, param = ScanBamParam(what = read_name_sam))))

  unique_input_reads <- unique(input_reads)
  n_unique_input_reads <- length(unique_input_reads)

  chosen_unique_output_reads_ix <- sample(1:n_unique_input_reads, size = n_unique_output_reads, replace = FALSE)
  chosen_unique_output_reads <- unique_input_reads[chosen_unique_output_reads_ix]
  filter <- FilterRules(list(filter = function(x) x$qname %in% chosen_unique_output_reads))
  if (is_transcriptome) {
    filterBam(
    file = input_name,
    destination = output_name,
    filter = filter,
    index = character(0),
    indexDestination = FALSE,
    param = ScanBamParam(what = scanBamWhat())
  )
  }
  else {
    filterBam(
    file = input_name,
    destination = output_name,
    filter = filter,
    indexDestination = TRUE,
    param = ScanBamParam(what = scanBamWhat())
  )
  }

  return(0)
}

suppressPackageStartupMessages(library(Rsamtools))

status <- subsample_bam(
  input_name = snakemake@input[["input_name"]],
  output_name = snakemake@output[["output_name"]],
  n_unique_output_reads = snakemake@params[["number_to_sample"]],
  seed = snakemake@params[["seed"]],
  read_name_sam = snakemake@params[["read_name_sam"]],
  is_transcriptome = snakemake@params[["is_transcriptome"]]
)

sessionInfo()

sink()
sink()

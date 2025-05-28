make_table_S01 <- function(output_path) {
  config <- rjson::fromJSON(file = "config/config.json")
  params <- rjson::fromJSON(file = "config/params.json")

  novel_counts <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method_novel, "_novel/gencode/transcript_counts_formatted.tsv")
  )[, -7]
  transcriptome_coding <- import(
    "results/annotate/run_orfanage/transcriptome.gtf"
  )

  novel_counts[, -(1:2)] <- edgeR::cpm(novel_counts[, -(1:2)])
  novel_counts$coding <- FALSE
  coding_transcripts <- unique(transcriptome_coding[transcriptome_coding$type == "CDS"]$transcript_id)

  novel_counts[which(novel_counts$transcript_id %in% coding_transcripts), ]$coding <- TRUE
  sqanti_categories <- vroom::vroom("results/annotate/run_sqanti/kinnex_wtc_11/kinnex_wtc_11_classification.txt") %>%
    dplyr::select(isoform, structural_category, subcategory) %>%
    rename(isoform = "transcript_id")

  novel_counts %>%
    left_join(sqanti_categories) %>%
    dplyr::select(
      gene_id,
      transcript_id,
      coding,
      structural_category,
      subcategory,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day1-rep1`,
      `day3-rep1`,
      `day3-rep2`,
      `day4-rep1`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    ) %>%
    write_tsv(output_path)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(rjson)
  library(dplyr)
  library(edgeR)
  library(vroom)
  library(readr)
  library(rtracklayer)
})

sessionInfo()
status <- make_table_S01(
  output_path = snakemake@output[[1]]
)

sink()
sink()

make_table_S03 <- function(output_path) {
  config <- rjson::fromJSON(file = "config/config.json")
  params <- rjson::fromJSON(file = "config/params.json")

  novel_counts <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method_novel, "_novel/gencode/transcript_counts_formatted.tsv")
  )

  y <- DGEList(
    counts = novel_counts[, -c(1:2, 7)],
    genes = data.frame(tx_id = novel_counts$transcript_id, gene_id = novel_counts$gene_id)
  )
  y <- calcNormFactors(y)

  keep <- rowSums(edgeR::cpm(y) > config$cpm_threshold) >= config$cpm_sample_threshold
  y <- y[keep, ]
  cps <- edgeR::cpm(y)

  rep <- c(
    "rep1",
    "rep2",
    "rep3",
    "rep1",
    "rep1",
    "rep2",
    "rep1",
    "rep1",
    "rep2",
    "rep3"
  )
  day <- c(
    rep(0, 3),
    1,
    rep(3, 2),
    4,
    rep(5, 3)
  )
  md <- data.frame(
    rep = rep,
    day = day
  )
  mm <- model.matrix(~ rep + ns(day, df = 2), data = md)
  y <- estimateDisp(y, design = mm)
  y <- glmFit(y, design = mm, prior.count = 1)
  ds <- diffSpliceDGE(y, geneid = "gene_id", contrast = c(0,0,0,1,1))
  topSpliceDGE(ds, number = Inf) %>% write_tsv(output_path)
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
  library(limma)
  library(biomaRt)
  library(splines)
})

sessionInfo()
status <- make_table_S03(
  output_path = snakemake@output[[1]]
)

sink()
sink()

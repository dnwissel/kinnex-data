make_table_S02 <- function(output_path) {
  config <- rjson::fromJSON(file = "config/config.json")
  params <- rjson::fromJSON(file = "config/params.json")
  
  novel_counts <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method_novel, "_novel/gencode/gene_counts_formatted.tsv")
  )

  mart <- useEnsembl("ensembl", "hsapiens_gene_ensembl")
  z <- getBM(c("ensembl_gene_id", "hgnc_symbol"), "ensembl_gene_id", sapply(strsplit(novel_counts$gene_id, "\\."), function(x) x[[1]]), mart)

  z$hgnc_symbol[which(z$hgnc_symbol == "")] <- z$ensembl_gene_id[which(z$hgnc_symbol == "")]

  gene_name <- data.frame(
    gene_id = sapply(strsplit(novel_counts$gene_id, "\\."), function(x) x[[1]])
  ) %>%
    left_join(data.frame(
      gene_id = z$ensembl_gene_id,
      gene_name = z$hgnc_symbol
    ), multiple = "first")


  y <- DGEList(
    counts = novel_counts[, c(2:4, 10:12)],
    genes = data.frame(gene_id = gene_name$gene_id, gene_name = gene_name$gene_name)
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
    "rep2",
    "rep3"
  )
  day5 <- c(
    rep(0, 3),
    rep(1, 3)
  )
  mm <- model.matrix(~ day5 + rep)
  dge <- edgeR::estimateDisp(y, mm)
  fit <- edgeR::glmQLFit(dge, mm)
  mc <- limma::makeContrasts(day5, levels = colnames(fit$coefficients))
  qlf <- edgeR::glmQLFTest(fit, contrast = mc)
  tt <- edgeR::topTags(qlf, n = Inf)

  data.frame(tt) %>% write_tsv(output_path)
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
})

sessionInfo()
status <- make_table_S02(
  output_path = snakemake@output[[1]]
)

sink()
sink()

plot_figure_02 <- function(output_path, depth) {
  config <- rjson::fromJSON(file = "config/config.json")
  cols <- carto_pal(config$n_colors, "Safe")
  replicate_cols <- cols[-config$exclude_color]
  replicate_cols <- c(replicate_cols[1], replicate_cols[5], replicate_cols[3], replicate_cols[4], replicate_cols[2])

  replicate_cols <- c(
    # Day 0-1
    adjustcolor(replicate_cols[1], alpha.f = 1),
    adjustcolor(replicate_cols[1], alpha.f = 0.75),
    adjustcolor(replicate_cols[1], alpha.f = 0.5),
    # Day 1-1
    adjustcolor(replicate_cols[2], alpha.f = 1),
    # Day 3-1
    adjustcolor(replicate_cols[3], alpha.f = 1),
    adjustcolor(replicate_cols[3], alpha.f = 0.75),
    # Day 4-1
    adjustcolor(replicate_cols[4], alpha.f = 1),
    # Day 5-1
    adjustcolor(replicate_cols[5], alpha.f = 1),
    adjustcolor(replicate_cols[5], alpha.f = 0.75),
    adjustcolor(replicate_cols[5], alpha.f = 0.5)
  )


  coldata_ill <- data.frame(
    files = c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3/quant.sf")
    ),
    names = c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3",
      "day5-rep1",
      "day5-rep2",
      "day5-rep3"
    )
  )

  se <- tximeta(coldata_ill)
  illumina_se <- fishpond::computeInfRV(se)

  coldata_kinnex <- data.frame(
    files = c(
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day0-rep1/day0-rep1.quant"),
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day0-rep2/day0-rep2.quant"),
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day0-rep3/day0-rep3.quant"),
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day5-rep1/day5-rep1.quant"),
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day5-rep2/day5-rep2.quant"),
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day5-rep3/day5-rep3.quant")
    ),
    names = c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3",
      "day5-rep1",
      "day5-rep2",
      "day5-rep3"
    )
  )

  se <- tximeta(coldata_kinnex, type = "oarfish")

  kinnex_se <- fishpond::computeInfRV(se)
  novel_transcriptome <- import("results/discover/fix_bambu_gene_ids/transcriptome.gtf")
  gene_mapping <- data.frame(
    novel_transcriptome
  ) %>%
    filter(type %in% c("transcript", "gene")) %>%
    dplyr::select(transcript_id, gene_id) %>%
    group_by(gene_id) %>%
    mutate(n = n())

  plt_frame <- data.frame(rowData(illumina_se)) %>%
    rownames_to_column() %>%
    rename(meanInfRV = "illumina") %>%
    left_join(
      (data.frame(rowData(kinnex_se)) %>% rownames_to_column())
    ) %>%
    replace(is.na(.), 0) %>%
    rename(rowname = "transcript_id") %>%
    left_join(gene_mapping) %>%
    rename(meanInfRV = "kinnex")

  plt_frame$n <- ifelse(
    plt_frame$n < 5,
    as.character(plt_frame$n),
    "5+"
  )

  plt_frame$n <- factor(
    plt_frame$n,
    levels = c(as.character(1:4), "5+")
  )
  colnames(plt_frame)[2:3] <- c("Illumina", "Kinnex")
  set.seed(1)
  inf_var <- plt_frame %>%
    pivot_longer(cols = c("Illumina", "Kinnex"), values_to = "mean_inf_rv", names_to = "tech") %>%
    mutate(
      tech = recode(
        tech,
        `Illumina` = "Illumina",
        `Kinnex` = "Kinnex"
      )
    ) %>%
    ggplot(aes(x = n, y = mean_inf_rv)) +
    geom_violin(scale = "width") +
    scale_y_log10() +
    facet_wrap(~tech) +
    ggpubfigs::theme_big_simple() +
    geom_hline(yintercept = 1, lty = 2, linewidth = 1, color = "red") +
    labs(x = "Transcripts per gene (GENCODE)", y = "Inferential variability")


  all_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("salmon_illumina", "oarfish"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/gencode/", rep, "/30000000.0/transcript_counts_formatted.tsv")
        ) %>%
          arrange(desc(transcript_id))
      })
    }
  ))


  flips <- data.frame(
    tech = rep(c("Illumina", "Kinnex"), each = 12),
    transcript_id = rep(
      rep(
        c(
          "ENST00000400585",
          "ENST00000262608",
          "ENST00000612582",
          "ENST00000355219"
        ),
        each = 3
      ), 2
    ),
    cpm = c(
      edgeR::cpm(all_counts[[1]][, -(1:2)])[which((all_counts[[1]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000099954.19"]$transcript_id)), ][3, 1:3],
      edgeR::cpm(all_counts[[1]][, -(1:2)])[which((all_counts[[1]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000099954.19"]$transcript_id)), ][6, 1:3],
      edgeR::cpm(all_counts[[1]][, -(1:2)])[which((all_counts[[1]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000099954.19"]$transcript_id)), ][1, 1:3],
      edgeR::cpm(all_counts[[1]][, -(1:2)])[which((all_counts[[1]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000099954.19"]$transcript_id)), ][4, 1:3],
      edgeR::cpm(all_counts[[2]][, -(1:2)])[which((all_counts[[2]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000099954.19"]$transcript_id)), ][3, 1:3],
      edgeR::cpm(all_counts[[2]][, -(1:2)])[which((all_counts[[2]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000099954.19"]$transcript_id)), ][6, 1:3],
      edgeR::cpm(all_counts[[2]][, -(1:2)])[which((all_counts[[2]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000099954.19"]$transcript_id)), ][1, 1:3],
      edgeR::cpm(all_counts[[2]][, -(1:2)])[which((all_counts[[2]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000099954.19"]$transcript_id)), ][4, 1:3]
    ),
    rep = rep(
      rep(c("Day0-1", "Day0-2", "Day0-3"), 4),
      2
    )
  ) %>% ggplot(aes(x = transcript_id, fill = tech, y = cpm + 1)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~rep, nrow = 1) +
    theme_big_simple() +
    scale_fill_manual(
      values = c("#ffb441", "#D8178C")
    ) +
    scale_y_log10() +
    labs(y = "CPM", x = "", fill = "") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_blank())


  flip_track <- (ggdraw() +
    draw_image(magick::image_read_pdf(config$flip_track_names,
      density = config$pdf_density
    ))
  )


  all_counts <- Reduce(c, lapply(
    c("oarfish", "salmon_illumina"), function(method) {
      lapply(
        c("5000000.0", "10000000.0", "20000000.0", "30000000.0"), function(depth) {
          vroom::vroom(
            paste0("results/format_quantify_subsampled/", method, "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
          )
        }
      )
    }
  ))

  joint_isoforms <- Reduce(intersect, lapply(all_counts, function(counts) {
    which(apply(edgeR::cpm(counts[, -(1:2)]), 1, function(x) sum(x > 1.0) >= 3))
  }))

  all_counts_filtered <- lapply(
    all_counts, function(x) x[joint_isoforms, ]
  )

  dte_calls <- sapply(
    all_counts,
    function(x) {
      dge <- edgeR::DGEList(counts = x[, -(1:2)], genes = x[, 2], group = c(0, 0, 0, 1, 1, 1))


      grp <- sapply(colnames(dge), function(x) strsplit(x, "-")[[1]][1])
      rep <- sapply(colnames(dge), function(x) strsplit(x, "-")[[1]][2])
      mm <- model.matrix(~ grp + rep)
      keep <- edgeR::filterByExpr(dge, mm)
      dge <- dge[keep, ]

      dge <- edgeR::estimateDisp(dge, mm)

      fit <- edgeR::glmQLFit(dge, mm)
      mc <- limma::makeContrasts(grpday5, levels = colnames(fit$coefficients))

      qlf <- edgeR::glmQLFTest(fit, contrast = mc)
      tt <- edgeR::topTags(qlf, n = Inf)
      return(tt$table$transcript_id[which(tt$table$FDR < 0.01)])
    }
  )

  all_calls <- unique(unlist(dte_calls))

  plt_frame <- data.frame(sapply(
    dte_calls,
    function(x) {
      as.numeric(all_calls %in% x)
    }
  ))

  colnames(plt_frame) <- paste0(rep(c(
    "Kinnex ",
    "Illumina "
  ), each = 4), rep(c("5 M", "10 M", "20 M", "30 M"), 2))



  dte_gencode <- upset(
    data = plt_frame,
    intersect = rev(c(
      "Kinnex 30 M",
      "Kinnex 20 M",
      "Kinnex 10 M",
      "Kinnex 5 M",
      "Illumina 30 M",
      "Illumina 20 M",
      "Illumina 10 M",
      "Illumina 5 M"
    )),
    width_ratio = 0.1, min_size = 25, wrap = TRUE, set_sizes = FALSE,
    sort_sets = FALSE,
    themes = upset_default_themes(text = element_text(size = 16)),
    base_annotations = list("Intersection size" = intersection_size(counts = FALSE)),
  ) + labs(title = "    DTE (GENCODE) (no correction)") + theme(plot.title = element_text(size = 20))

  all_counts <- Reduce(c, lapply(
    c("oarfish", "run_salmon_illumina_bootstrap"), function(method) {
      lapply(
        c("5000000.0", "10000000.0", "20000000.0", "30000000.0"), function(depth) {
          if (method == "run_salmon_illumina_bootstrap") {
            s <- edgeR::catchSalmon(
              c(
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
              )
            )
            sort_ix <- sort(rownames(s$annotation), index.return = TRUE)$ix
            DGEList(counts = (s$counts / s$annotation$Overdispersion)[sort_ix, ], genes = data.frame(transcript_id = rownames(s$annotation)[sort_ix]))
          } else {
            s <- vroom::vroom(
              paste0("results/format_quantify_subsampled/", method, "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
            )

            annot <- data.frame(
              transcript_id = s$transcript_id,
              gene_id = s$gene_id
            )

            rownames(annot) <- annot$transcript_id
            sort_ix <- sort(rownames(annot), index.return = TRUE)$ix
            DGEList(
              counts = s[sort_ix, -(1:2)],
              genes = annot[sort_ix, ]
            )
          }
        }
      )
    }
  ))

  joint_isoforms <- Reduce(intersect, lapply(all_counts, function(counts) {
    which(apply(edgeR::cpm(counts$counts), 1, function(x) sum(x > 1.0) >= 3))
  }))

  all_counts_filtered <- lapply(
    all_counts, function(x) x[joint_isoforms, ]
  )

  dte_calls <- sapply(
    all_counts,
    function(x) {
      dge <- x
      grp <- factor(c(rep(0, 3), rep(1, 3)))
      rep <- factor(rep(1:3, 2))
      mm <- model.matrix(~ grp + rep)
      keep <- filterByExpr(dge, mm)
      dge <- dge[keep, ]
      dge <- edgeR::estimateDisp(dge, mm)

      fit <- edgeR::glmQLFit(dge, mm)
      mc <- limma::makeContrasts(grp1, levels = colnames(fit$coefficients))

      qlf <- edgeR::glmQLFTest(fit, contrast = mc)
      tt <- edgeR::topTags(qlf, n = Inf)
      return(tt$table$transcript_id[which(tt$table$FDR < 0.01)])
    }
  )

  all_calls <- unique(unlist(dte_calls))

  plt_frame <- data.frame(sapply(
    dte_calls,
    function(x) {
      as.numeric(all_calls %in% x)
    }
  ))

  colnames(plt_frame) <- paste0(rep(c(
    "Kinnex ",
    "Illumina "
  ), each = 4), rep(c("5 M", "10 M", "20 M", "30 M"), 2))



  dte_gencode_corrected <- upset(
    data = plt_frame,
    intersect = rev(c(
      "Kinnex 30 M",
      "Kinnex 20 M",
      "Kinnex 10 M",
      "Kinnex 5 M",
      "Illumina 30 M",
      "Illumina 20 M",
      "Illumina 10 M",
      "Illumina 5 M"
    )),
    width_ratio = 0.1, min_size = 100, wrap = TRUE, set_sizes = FALSE,
    sort_sets = FALSE,
    themes = upset_default_themes(text = element_text(size = 16)),
    base_annotations = list("Intersection size" = intersection_size(counts = FALSE)),
  ) + labs(title = "    DTE (GENCODE) (with correction)") + theme(plot.title = element_text(size = 20))

  transcriptome <- import("results/discover/fix_bambu_gene_ids/transcriptome.gtf")

  length_frame <- data.frame(transcriptome) %>%
    filter(type == "exon") %>%
    mutate(exonic_length = end - start) %>%
    group_by(transcript_id) %>%
    summarise(exonic_length = sum(exonic_length))

  kinnex_frame <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", config$chosen_method, "/gencode/1/30000000.0/transcript_counts_formatted.tsv")
  ) %>% arrange(desc(transcript_id))

  s <- edgeR::catchSalmon(
    c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
    )
  )

  illumina_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day5-rep1` = s$counts[, 4],
    `day5-rep2` = s$counts[, 5],
    `day5-rep3` = s$counts[, 6],
    check.names = FALSE
  ) %>%
    arrange(desc(transcript_id)) %>%
    left_join(
      (
        kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
      )
    ) %>%
    dplyr::select(
      gene_id, transcript_id,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    )


  kinnex_frame_gene <- kinnex_frame %>%
    group_by(gene_id) %>%
    summarise(
      `day0-rep1` = sum(`day0-rep1`),
      `day0-rep2` = sum(`day0-rep2`),
      `day0-rep3` = sum(`day0-rep3`),
      `day5-rep1` = sum(`day5-rep1`),
      `day5-rep2` = sum(`day5-rep2`),
      `day5-rep3` = sum(`day5-rep3`)
    )

  illumina_frame_gene <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", "salmon_illumina", "/gencode/1/30000000.0/transcript_counts_formatted.tsv")
  ) %>%
    arrange(desc(transcript_id)) %>%
    group_by(gene_id) %>%
    summarise(
      `day0-rep1` = sum(`day0-rep1`),
      `day0-rep2` = sum(`day0-rep2`),
      `day0-rep3` = sum(`day0-rep3`),
      `day5-rep1` = sum(`day5-rep1`),
      `day5-rep2` = sum(`day5-rep2`),
      `day5-rep3` = sum(`day5-rep3`)
    )

  n_kinnex_transcripts <- apply(edgeR::cpm(kinnex_frame[, -(1:2)]), 2, function(x) rep(1, sum(x > config$cpm_threshold)))
  n_illumina_transcripts <- apply(edgeR::cpm(illumina_frame[, -(1:2)]), 2, function(x) rep(1, sum(x > config$cpm_threshold)))

  n_kinnex_genes <- apply(edgeR::cpm(kinnex_frame_gene[, -1]), 2, function(x) rep(1, sum(x > config$cpm_threshold)))
  n_illumina_genes <- apply(edgeR::cpm(illumina_frame_gene[, -1]), 2, function(x) rep(1, sum(x > config$cpm_threshold)))

  plt_frame <- data.frame(
    n = c(
      unlist(n_kinnex_transcripts),
      unlist(n_illumina_transcripts),
      unlist(n_kinnex_genes),
      unlist(n_illumina_genes)
    ),
    tech = c(
      rep("Kinnex", length(unlist(n_kinnex_transcripts))),
      rep("Illumina", length(unlist(n_illumina_transcripts))),
      rep("Kinnex", length(unlist(n_kinnex_genes))),
      rep("Illumina", length(unlist(n_illumina_genes)))
    ),
    type = c(
      rep("Transcript", length(unlist(n_kinnex_transcripts))),
      rep("Transcript", length(unlist(n_illumina_transcripts))),
      rep("Gene", length(unlist(n_kinnex_genes))),
      rep("Gene", length(unlist(n_illumina_genes)))
    ),
    day = c(
      rep(c("Day0-1", "Day0-2", "Day0-3", "Day5-1", "Day5-2", "Day5-3"), sapply(n_kinnex_transcripts, length)),
      rep(c("Day0-1", "Day0-2", "Day0-3", "Day5-1", "Day5-2", "Day5-3"), sapply(n_illumina_transcripts, length)),
      rep(c("Day0-1", "Day0-2", "Day0-3", "Day5-1", "Day5-2", "Day5-3"), sapply(n_kinnex_genes, length)),
      rep(c("Day0-1", "Day0-2", "Day0-3", "Day5-1", "Day5-2", "Day5-3"), sapply(n_illumina_genes, length))
    )
  )

  transcripts_thirty_downsampled <- plt_frame %>%
    filter(type == "Transcript") %>%
    ggplot(aes(x = day)) +
    geom_bar(position = "stack") +
    facet_wrap(~tech, nrow = 1) +
    theme_big_simple() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = "", x = "Sample", y = "Transcripts") +
    scale_fill_manual(values = rev(config$known_novel_cols))

  genes_thirty_downsampled <- plt_frame %>%
    filter(type == "Gene") %>%
    ggplot(aes(x = day)) +
    geom_bar(position = "stack") +
    facet_wrap(~tech, nrow = 1) +
    theme_big_simple() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = "", x = "Sample", y = "Genes") +
    scale_fill_manual(values = rev(config$known_novel_cols))

  identified_day_zero <- lapply(c("30000000.0", "20000000.0", "10000000.0", "5000000.0"), function(depth) {
    lapply(c("salmon_illumina", config$chosen_method), function(method) {
      gene_frame <- vroom::vroom(
        paste0("results/format_quantify_subsampled/", method, "/gencode/1/", depth, "/gene_counts_formatted.tsv")
      ) %>% arrange(desc(gene_id))

      gene_frame$gene_id[which(apply(edgeR::cpm(gene_frame[, -1])[, 1:3], 1, function(x) (sum(x >= config$cpm_threshold) >= 3)))]
    })
  })

  all_genes_identified <- unique(unlist(identified_day_zero))

  identified_day_zero <- list(
    identified_day_zero[[1]][[2]],
    identified_day_zero[[2]][[2]],
    identified_day_zero[[3]][[2]],
    identified_day_zero[[4]][[2]],
    identified_day_zero[[1]][[1]],
    identified_day_zero[[2]][[1]],
    identified_day_zero[[3]][[1]],
    identified_day_zero[[4]][[1]]
  )

  plt_frame <- data.frame(sapply(
    identified_day_zero,
    function(x) {
      as.numeric(all_genes_identified %in% x)
    }
  ))

  colnames(plt_frame) <- paste0(rep(c(
    "Kinnex ",
    "Illumina "
  ), each = 4), rep(rev(c("5 M", "10 M", "20 M", "30 M")), 2))



  identified_gene <- upset(
    data = plt_frame,
    intersect = rev(c(
      "Kinnex 30 M",
      "Kinnex 20 M",
      "Kinnex 10 M",
      "Kinnex 5 M",
      "Illumina 30 M",
      "Illumina 20 M",
      "Illumina 10 M",
      "Illumina 5 M"
    )),
    width_ratio = 0.1, min_size = 50, wrap = TRUE, set_sizes = FALSE,
    sort_sets = FALSE,
    themes = upset_default_themes(text = element_text(size = 16)),
    base_annotations = list("Intersection size" = intersection_size(counts = FALSE)),
  ) + labs(title = "    Identified genes (Day 0)") + theme(plot.title = element_text(size = 20))


  identified_day_zero <- lapply(c("30000000.0", "20000000.0", "10000000.0", "5000000.0"), function(depth) {
    lapply(c("salmon_illumina", "oarfish"), function(method) {
      if (method == "oarfish") {
        transcript_frame <- vroom::vroom(
          paste0("results/format_quantify_subsampled/", "oarfish", "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
        ) %>% arrange(desc(transcript_id))
      } else {
        s <- edgeR::catchSalmon(
          c(
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
          )
        )

        transcript_frame <- data.frame(
          transcript_id = rownames(s$counts),
          `day0-rep1` = s$counts[, 1],
          `day0-rep2` = s$counts[, 2],
          `day0-rep3` = s$counts[, 3],
          `day5-rep1` = s$counts[, 4],
          `day5-rep2` = s$counts[, 5],
          `day5-rep3` = s$counts[, 6],
          check.names = FALSE
        ) %>%
          arrange(desc(transcript_id)) %>%
          left_join(
            (
              kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
            )
          ) %>%
          dplyr::select(
            gene_id, transcript_id,
            `day0-rep1`,
            `day0-rep2`,
            `day0-rep3`,
            `day5-rep1`,
            `day5-rep2`,
            `day5-rep3`
          )
      }

      transcript_frame$transcript_id[which(apply(edgeR::cpm(transcript_frame[, -(1:2)])[, 1:3], 1, function(x) (sum(x >= config$cpm_threshold) >= 3)))]
    })
  })

  all_transcripts_identified <- unique(unlist(identified_day_zero))

  identified_day_zero <- list(
    identified_day_zero[[1]][[2]],
    identified_day_zero[[2]][[2]],
    identified_day_zero[[3]][[2]],
    identified_day_zero[[4]][[2]],
    identified_day_zero[[1]][[1]],
    identified_day_zero[[2]][[1]],
    identified_day_zero[[3]][[1]],
    identified_day_zero[[4]][[1]]
  )

  plt_frame <- data.frame(sapply(
    identified_day_zero,
    function(x) {
      as.numeric(all_transcripts_identified %in% x)
    }
  ))

  colnames(plt_frame) <- paste0(rep(c(
    "Kinnex ",
    "Illumina "
  ), each = 4), rep(rev(c("5 M", "10 M", "20 M", "30 M")), 2))



  identified_transcripts <- upset(
    data = plt_frame,
    intersect = rev(c(
      "Kinnex 30 M",
      "Kinnex 20 M",
      "Kinnex 10 M",
      "Kinnex 5 M",
      "Illumina 30 M",
      "Illumina 20 M",
      "Illumina 10 M",
      "Illumina 5 M"
    )),
    width_ratio = 0.1, min_size = 150, wrap = TRUE, set_sizes = FALSE,
    sort_sets = FALSE,
    themes = upset_default_themes(text = element_text(size = 16)),
    base_annotations = list("Intersection size" = intersection_size(counts = FALSE)),
  ) + labs(title = "    Identified transcripts (Day 0)") + theme(plot.title = element_text(size = 20))




  all_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("salmon_illumina", "oarfish"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/gencode/", rep, "/30000000.0/transcript_counts_formatted.tsv")
        ) %>%
          arrange(desc(transcript_id))
      })
    }
  ))

  division <- data.frame(
    tech = c(rep("Illumina", 12), rep("Kinnex", 12)),
    transcript_id = c(
      rep("ENST00000647818", 3),
      rep("ENST00000603108", 3),
      rep("ENST00000648964", 3),
      rep("ENST00000263440", 3),
      rep("ENST00000647818", 3),
      rep("ENST00000603108", 3),
      rep("ENST00000648964", 3),
      rep("ENST00000263440", 3)
    ),
    cpm = c(
      edgeR::cpm(all_counts[[1]][, -(1:2)])[which((all_counts[[1]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000088035.18"]$transcript_id)), ][6, 1:3],
      edgeR::cpm(all_counts[[1]][, -(1:2)])[which((all_counts[[1]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000088035.18"]$transcript_id)), ][7, 1:3],
      edgeR::cpm(all_counts[[1]][, -(1:2)])[which((all_counts[[1]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000088035.18"]$transcript_id)), ][5, 1:3],
      edgeR::cpm(all_counts[[1]][, -(1:2)])[which((all_counts[[1]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000088035.18"]$transcript_id)), ][12, 1:3],
      edgeR::cpm(all_counts[[2]][, -(1:2)])[which((all_counts[[2]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000088035.18"]$transcript_id)), ][6, 1:3],
      edgeR::cpm(all_counts[[2]][, -(1:2)])[which((all_counts[[2]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000088035.18"]$transcript_id)), ][7, 1:3],
      edgeR::cpm(all_counts[[2]][, -(1:2)])[which((all_counts[[2]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000088035.18"]$transcript_id)), ][5, 1:3],
      edgeR::cpm(all_counts[[2]][, -(1:2)])[which((all_counts[[2]])$transcript_id %in% unique(novel_transcriptome[novel_transcriptome$gene_id == "ENSG00000088035.18"]$transcript_id)), ][12, 1:3]
    ),
    rep = rep(c("Day0-1", "Day0-2", "Day0-3"), 4 * 2)
  ) %>% ggplot(aes(x = transcript_id, fill = tech, y = cpm + 1)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
    facet_wrap(~rep, nrow = 1) +
    theme_big_simple() +
    scale_fill_manual(
      values = c("#ffb441", "#D8178C")
    ) +
    scale_y_log10() +
    labs(y = "CPM", x = "", fill = "") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_blank())

  identified_day_zero <- lapply(c("30000000.0"), function(depth) {
    lapply(c("salmon_illumina", config$chosen_method), function(method) {
      if (method == config$chosen_method) {
        transcript_frame <- vroom::vroom(
          paste0("results/format_quantify_subsampled/", config$chosen_method, "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
        ) %>% arrange(desc(transcript_id))
      } else {
        s <- edgeR::catchSalmon(
          c(
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
            paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
          )
        )

        transcript_frame <- data.frame(
          transcript_id = rownames(s$counts),
          `day0-rep1` = s$counts[, 1],
          `day0-rep2` = s$counts[, 2],
          `day0-rep3` = s$counts[, 3],
          `day5-rep1` = s$counts[, 4],
          `day5-rep2` = s$counts[, 5],
          `day5-rep3` = s$counts[, 6],
          check.names = FALSE
        ) %>%
          arrange(desc(transcript_id)) %>%
          left_join(
            (
              kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
            )
          ) %>%
          dplyr::select(
            gene_id, transcript_id,
            `day0-rep1`,
            `day0-rep2`,
            `day0-rep3`,
            `day5-rep1`,
            `day5-rep2`,
            `day5-rep3`
          )
      }

      transcript_frame$transcript_id[which(apply(edgeR::cpm(transcript_frame[, -(1:2)])[, 1:3], 1, function(x) (sum(x >= 1) >= 3)))]
    })
  })[[1]]



  s <- edgeR::catchSalmon(
    c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
    )
  )

  illumina_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day5-rep1` = s$counts[, 4],
    `day5-rep2` = s$counts[, 5],
    `day5-rep3` = s$counts[, 6],
    check.names = FALSE
  ) %>%
    arrange(desc(transcript_id)) %>%
    left_join(
      (
        kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
      )
    ) %>%
    dplyr::select(
      gene_id, transcript_id,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    )

  plt_frame <- data.frame(
    transcript_id = c(
      intersect(identified_day_zero[[1]], identified_day_zero[[2]]),
      setdiff(identified_day_zero[[1]], identified_day_zero[[2]]),
      setdiff(identified_day_zero[[2]], identified_day_zero[[1]])
    ),
    group = c(
      rep("Both", length(intersect(identified_day_zero[[1]], identified_day_zero[[2]]))),
      rep("Illumina", length(setdiff(identified_day_zero[[1]], identified_day_zero[[2]]))),
      rep("Kinnex", length(setdiff(identified_day_zero[[2]], identified_day_zero[[1]])))
    )
  ) %>% left_join(
    data.frame(
      transcript_id = illumina_frame$transcript_id,
      rep_var = apply(edgeR::cpm(illumina_frame[, 3:5]), 1, function(x) var(x))
    )
  )

  transcriptome_sequences <- readDNAStringSet("results/prepare/extract_transcriptomes/gencode_transcriptome.fa")
  names(transcriptome_sequences) <- sapply(strsplit(names(transcriptome_sequences), "\\ "), function(x) x[[1]])

  plt_frame <- plt_frame %>% left_join(
    data.frame(
      transcript_id = names(transcriptome_sequences),
      gc_content = letterFrequency(transcriptome_sequences, letters = "GC", as.prob = TRUE)[, 1],
      length = nchar(transcriptome_sequences)
    )
  )

  coldata_ill <- data.frame(
    files = c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3/quant.sf")
    ),
    names = c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3"
    )
  )

  se <- tximeta(coldata_ill)
  illumina_se <- fishpond::computeInfRV(se)

  plt_frame <- plt_frame %>% left_join(
    data.frame(transcript_id = names(rowData(illumina_se)[, 1]), inf_rv = unname(rowData(illumina_se)[, 1]))
  )

  transcript_reasoning <- plt_frame %>%
    ggplot(aes(x = group, y = length)) +
    geom_violin() +
    geom_hline(yintercept = 1250, lty = 2, color = "red", linewidth = 1) +
    scale_y_log10() +
    theme_big_simple() +
    labs(x = "Detected by", y = "Transcript length", fill = "") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  fourth_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        (ggdraw() +
          draw_image(magick::image_read_pdf(config$division_track_adapt,
            density = config$pdf_density
          ))
        ),
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 1,
    labels = c("H"),
    label_size = 36
  )

  illumina_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day5-rep1` = s$counts[, 4],
    `day5-rep2` = s$counts[, 5],
    `day5-rep3` = s$counts[, 6],
    check.names = FALSE
  ) %>%
    arrange(desc(transcript_id)) %>%
    left_join(
      (
        kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
      )
    ) %>%
    dplyr::select(
      gene_id, transcript_id,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    )

  kinnex_frame <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", "bambu", "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
  ) %>% arrange(desc(transcript_id))

  s <- edgeR::catchSalmon(
    c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
    )
  )

  illumina_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day5-rep1` = s$counts[, 4],
    `day5-rep2` = s$counts[, 5],
    `day5-rep3` = s$counts[, 6],
    check.names = FALSE
  ) %>%
    arrange(desc(transcript_id)) %>%
    left_join(
      (
        kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
      )
    ) %>%
    dplyr::select(
      gene_id, transcript_id,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    )

  kinnex_frame <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", "oarfish", "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
  ) %>% arrange(desc(transcript_id))

  plt_frame$reasoning <- "Unclear"

  plt_frame$reasoning <- ifelse(
    plt_frame$group == "Illumina" & plt_frame$length < 1250,
    "Length",
    "Unclear"
  )

  coldata_ill <- data.frame(
    files = c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3/quant.sf")
    ),
    names = c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3"
    )
  )

  se <- tximeta(coldata_ill)
  illumina_se <- fishpond::computeInfRV(se)

  coldata_kinnex <- data.frame(
    files = c(
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day0-rep1/day0-rep1.quant"),
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day0-rep2/day0-rep2.quant"),
      paste0("results/quantify_subsampled/run_oarfish_bootstrap/gencode/1/", depth, "/day0-rep3/day0-rep3.quant")
    ),
    names = c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3"
    )
  )

  se <- tximeta(coldata_kinnex, type = "oarfish")

  kinnex_se <- fishpond::computeInfRV(se)
  set.seed(config$seed)
  novel_transcriptome <- import("results/discover/fix_bambu_gene_ids/transcriptome.gtf")
  gene_mapping <- data.frame(
    novel_transcriptome
  ) %>%
    filter(type %in% c("transcript", "gene")) %>%
    dplyr::select(transcript_id, gene_id) %>%
    group_by(gene_id) %>%
    mutate(n = n())

  inf_var_plt_frame <- data.frame(rowData(illumina_se)) %>%
    rownames_to_column() %>%
    rename(meanInfRV = "illumina") %>%
    left_join(
      (data.frame(rowData(kinnex_se)) %>% rownames_to_column())
    ) %>%
    replace(is.na(.), 0) %>%
    rename(rowname = "transcript_id") %>%
    left_join(gene_mapping) %>%
    rename(meanInfRV = "kinnex")

  inf_var_plt_frame$n <- ifelse(
    inf_var_plt_frame$n < 7,
    as.character(inf_var_plt_frame$n),
    "7+"
  )

  inf_var_plt_frame$n <- factor(
    inf_var_plt_frame$n,
    levels = c(as.character(1:6), "7+")
  )
  colnames(inf_var_plt_frame)[2:3] <- c("Illumina", "Kinnex")

  kinnex_frame <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", config$chosen_method, "/gencode/1/30000000.0/transcript_counts_formatted.tsv")
  ) %>% arrange(desc(transcript_id))

  s <- edgeR::catchSalmon(
    c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
    )
  )

  illumina_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day5-rep1` = s$counts[, 4],
    `day5-rep2` = s$counts[, 5],
    `day5-rep3` = s$counts[, 6],
    check.names = FALSE
  ) %>%
    arrange(desc(transcript_id)) %>%
    left_join(
      (
        kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
      )
    ) %>%
    dplyr::select(
      gene_id, transcript_id,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    )


  illumina_multinomial <- illumina_frame[, 1:5] %>% arrange(transcript_id)
  kinnex_multinomial <- kinnex_frame[, 1:5] %>% arrange(transcript_id)

  ill_mask <- apply(edgeR::cpm(illumina_multinomial[, -(1:2)]), 1, function(x) all(x >= 1.0))
  kinnex_mask <- apply(edgeR::cpm(kinnex_multinomial[, -(1:2)]), 1, function(x) all(x >= 1.0))
  joint_mask <- ill_mask | kinnex_mask

  illumina_multinomial_filtered <- data.frame(illumina_multinomial[joint_mask, ], check.names = FALSE)
  kinnex_multinomial_filtered <- data.frame(kinnex_multinomial[joint_mask, ], check.names = FALSE)

  joint_genes <- intersect(illumina_multinomial_filtered$gene_id, kinnex_multinomial_filtered$gene_id)

  illumina_multinomial_filtered <- illumina_multinomial_filtered %>% filter(gene_id %in% joint_genes)
  kinnex_multinomial_filtered <- kinnex_multinomial_filtered %>% filter(gene_id %in% joint_genes)

  multi_genes <- names(which(table(illumina_multinomial_filtered$gene_id) > 1))

  illumina_cpm_rounded <- ((illumina_multinomial_filtered[, 3:5]))
  kinnex_cpm_rounded <- ((kinnex_multinomial_filtered[, 3:5]))

  dge <- DGEList(counts = cbind(as.matrix(illumina_cpm_rounded), as.matrix(kinnex_cpm_rounded)), genes = illumina_multinomial_filtered[, 1:2])

  grp <- factor(c(rep(0, 3), rep(1, 3)))
  rep <- factor(rep(1:3, 2))
  mm <- model.matrix(~ grp + rep)

  dge <- edgeR::estimateDisp(dge, mm)

  fit <- edgeR::glmQLFit(dge, mm)
  mc <- limma::makeContrasts(grp1, levels = colnames(fit$coefficients))
  ds <- diffSpliceDGE(fit, geneid = "gene_id", contrast = mc)

  tt <- topSpliceDGE(ds, number = Inf)

  tt_gene <- topSpliceDGE(ds, number = Inf, test = "gene")

  tt_exon <- topSpliceDGE(ds, number = Inf, test = "exon")

  p_values <- tt_exon %>% dplyr::select(transcript_id, FDR)



  inf_var_missed_transcripts <- inf_var_plt_frame %>%
    filter(transcript_id %in% (plt_frame %>% filter(group == "Kinnex") %>% pull(transcript_id)) & Illumina > 1) %>%
    pull(transcript_id)

  division_missed_transcripts <- p_values %>%
    filter(transcript_id %in% (plt_frame %>% filter(group == "Illumina") %>% pull(transcript_id)) & FDR <= 1e-05) %>%
    pull(transcript_id)

  plt_frame$reasoning[plt_frame$transcript_id %in% inf_var_missed_transcripts & (!plt_frame$transcript_id %in% (plt_frame %>% filter(reasoning == "Length") %>% pull(transcript_id)))] <- "Inf. var."

  plt_frame$reasoning[(plt_frame$transcript_id %in% division_missed_transcripts) & (!plt_frame$transcript_id %in% (plt_frame %>% filter(reasoning == "Length") %>% pull(transcript_id)))] <- "Division"

  s <- edgeR::catchSalmon(
    c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
    )
  )

  illumina_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day5-rep1` = s$counts[, 4],
    `day5-rep2` = s$counts[, 5],
    `day5-rep3` = s$counts[, 6],
    check.names = FALSE
  ) %>%
    arrange(desc(transcript_id)) %>%
    left_join(
      (
        kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
      )
    ) %>%
    dplyr::select(
      gene_id, transcript_id,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    )

  illumina_cpm <- edgeR::cpm(illumina_frame[, -(1:2)])

  kinnex_frame <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", "oarfish", "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
  ) %>% arrange(desc(transcript_id))

  kinnex_cpm <- edgeR::cpm(kinnex_frame[, -(1:2)])

  plt_frame <- plt_frame %>%
    filter(group != "Both") %>%
    left_join(rbind(
      data.frame(
        transcript_id = illumina_frame$transcript_id[match((plt_frame %>% filter(group == "Illumina") %>% pull(transcript_id)), illumina_frame$transcript_id)],
        mean_cpm = apply(illumina_cpm[, 1:3], 1, mean)[match((plt_frame %>% filter(group == "Illumina") %>% pull(transcript_id)), illumina_frame$transcript_id)]
      ),
      data.frame(
        transcript_id = kinnex_frame$transcript_id[match((plt_frame %>% filter(group == "Kinnex") %>% pull(transcript_id)), kinnex_frame$transcript_id)],
        mean_cpm = apply(kinnex_cpm[, 1:3], 1, mean)[match((plt_frame %>% filter(group == "Kinnex") %>% pull(transcript_id)), kinnex_frame$transcript_id)]
      )
    ), by = c("transcript_id" = "transcript_id"))

  plt_frame$cut_cpm <- cut(
    plt_frame$mean_cpm,
    c(0, 10, 100, 3225)
  )

  plt_frame$cut_cpm <- paste0("CPM: ", plt_frame$cut_cpm)

  plt_frame$cut_cpm <- ifelse(
    plt_frame$cut_cpm == "CPM: (0,10]",
    "CPM: (0,10]",
    ifelse(
      plt_frame$cut_cpm == "CPM: (10,100]",
      "CPM: (10,100]",
      "CPM > 100"
    )
  )
  plt_frame$cut_cpm <- factor(
    plt_frame$cut_cpm,
    levels = c(
      "CPM: (0,10]",
      "CPM: (10,100]",
      "CPM > 100"
    )
  )

  plt_frame$group <- ifelse(
    plt_frame$group == "Illumina",
    "Missed by Kinnex",
    "Missed by Illumina"
  )

  overall_counts <- plt_frame %>%
    ggplot(aes(x = reasoning)) +
    geom_bar() +
    facet_grid(cols = vars(cut_cpm), rows = vars(group), scales = "free_x") +
    ggpubfigs::theme_big_simple() +
    labs(x = "Reason for non-detection", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


  hm <- tximeta(coldata_ill)
  illumina_se <- fishpond::computeInfRV(hm)

  novel_transcriptome <- import("results/discover/fix_bambu_gene_ids/transcriptome.gtf")
  gene_mapping <- data.frame(
    novel_transcriptome
  ) %>%
    filter(type %in% c("transcript", "gene")) %>%
    dplyr::select(transcript_id, gene_id) %>%
    group_by(gene_id) %>%
    mutate(n = n())

  plt_frame <- data.frame(rowData(illumina_se)) %>%
    rownames_to_column() %>%
    rename(meanInfRV = "illumina") %>%
    replace(is.na(.), 0) %>%
    rename(rowname = "transcript_id") %>%
    left_join(gene_mapping)

  plt_frame$n <- ifelse(
    plt_frame$n < 7,
    as.character(plt_frame$n),
    "7+"
  )

  plt_frame$n <- factor(
    plt_frame$n,
    levels = c(as.character(1:6), "7+")
  )
  colnames(plt_frame)[2] <- c("Illumina")
  set.seed(42)

  non_variable_transcripts <- plt_frame %>%
    filter(Illumina <= 1) %>%
    pull(transcript_id)

  all_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("salmon_illumina", "oarfish"), function(method) {
        if (method == "salmon_illumina") {
          s <- edgeR::catchSalmon(
            c(
              paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
              paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
              paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
              paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
              paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
              paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
            )
          )

          data.frame(
            gene_id = NA,
            transcript_id = rownames(s$counts),
            `day0-rep1` = s$counts[, 1],
            `day0-rep2` = s$counts[, 2],
            `day0-rep3` = s$counts[, 3],
            `day5-rep1` = s$counts[, 4],
            `day5-rep2` = s$counts[, 5],
            `day5-rep3` = s$counts[, 6],
            check.names = FALSE
          ) %>%
            filter(transcript_id %in% non_variable_transcripts) %>%
            arrange(desc(transcript_id))
        } else {
          vroom::vroom(
            paste0("results/format_quantify_subsampled/", method, "/gencode/", rep, "/30000000.0/transcript_counts_formatted.tsv")
          ) %>%
            filter(transcript_id %in% non_variable_transcripts) %>%
            arrange(desc(transcript_id))
        }
      })
    }
  ))

  prepped_counts <- lapply(
    all_counts,
    function(x) {
      x <- DGEList(x[, -(1:2)], genes = x[, 1], group = c(rep("day0", 3), rep("day5", 3)))
      cpm <- edgeR::cpm(x, log = FALSE)
      return(cpm)
    }
  )

  mask <- intersect(
    which(apply(prepped_counts[[1]], 1, function(x) sum(x > 1.0) >= 6)),
    Reduce(union, lapply(prepped_counts[-1], function(counts) {
      which(apply(counts[, 1:3], 1, function(x) (sum(x > 1.0) >= 3)) | apply(counts[, 4:6], 1, function(x) (sum(x > 1.0) >= 3)))
    }))
  )

  mask_absolute <- intersect(
    which(apply(prepped_counts[[1]], 1, function(x) sum(x > 1.0) >= 6)),
    Reduce(union, lapply(prepped_counts[-1], function(counts) {
      which(apply(counts[, 1:3], 1, function(x) (sum(x > 1.0) >= 3)))
    }))
  )


  transcript_frame <- data.frame(transcript_id = all_counts[[1]]$transcript_id[mask]) %>% left_join(data.frame(novel_transcriptome) %>%
    mutate(total = end - start) %>%
    filter(type == "exon") %>% group_by(transcript_id) %>% summarise(length = sum(total)))

  prepped_counts <- lapply(
    all_counts,
    function(x) {
      x <- DGEList(x[mask, -(1:2)], genes = x[mask, 1], group = c(rep("day0", 3), rep("day5", 3)))
      cpm <- edgeR::cpm(x, log = TRUE)
      return(cpm)
    }
  )

  plt_frame <- data.frame(
    transcript_id = transcript_frame$transcript_id,
    length = transcript_frame$length,
    log_fc = unlist(lapply(prepped_counts[-1], function(counts) {
      c(
        counts[, 1] - counts[, 4],
        counts[, 2] - counts[, 5],
        counts[, 3] - counts[, 6]
      )
    })),
    illumina_fc = rep(
      c(
        prepped_counts[[1]][, 1] - prepped_counts[[1]][, 4],
        prepped_counts[[1]][, 2] - prepped_counts[[1]][, 5],
        prepped_counts[[1]][, 3] - prepped_counts[[1]][, 6]
      ),
      1
    ),
    rep = rep(rep(1:3, each = length(mask)), 1),
    method = rep(
      c(
        "Kinnex (Oarfish)"
      ),
      each = length(mask) * 3
    )
  )

  plt_frame$length_factor <- factor(cut(plt_frame$length, c(0, 750, 2000, 10000, 347559 + 1)))

  plt_frame$illumina_expr_factor <- factor(cut(plt_frame$illumina_fc, c(-9, -0.7412485, -0.1025366, 0.4798768, 8.3373789 + 1)))



  set.seed(config$seed)

  illumina_relative_transcript <- plt_frame %>%
    ggplot(aes(x = illumina_fc, y = log_fc, shape = as.factor(rep))) +
    geom_point(size = 4, alpha = .7) +
    labs(x = "Illumina (logFC)", y = "Kinnex (logFC)", shape = "Replicate") +
    ggpubfigs::theme_big_simple() +
    stat_cor(aes(label = ..r.label.., colour = NULL, group = 1), show.legend = FALSE, method = "pearson", size = 6, label.x.npc = c(0.025), label.y.npc = c(1)) +
    geom_density_2d(color = "lightblue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Relative quantification (transcript)")


  prepped_counts <- lapply(
    all_counts,
    function(x) {
      x <- DGEList(x[mask, -(1:2)], genes = x[mask, 1], group = c(rep("day0", 3), rep("day5", 3)))
      cpm <- edgeR::cpm(x, log = TRUE)
      return(cpm)
    }
  )

  plt_frame <- data.frame(
    transcript_id = all_counts[[1]]$transcript_id[mask],
    observed_abundance = unlist(lapply(prepped_counts[-1], function(counts) {
      c(
        counts[, 1],
        counts[, 2],
        counts[, 3],
        counts[, 4],
        counts[, 5],
        counts[, 6]
      )
    })),
    illumina_abundance = rep(
      c(
        prepped_counts[[1]][, 1],
        prepped_counts[[1]][, 2],
        prepped_counts[[1]][, 3],
        prepped_counts[[1]][, 4],
        prepped_counts[[1]][, 5],
        prepped_counts[[1]][, 6]
      ),
      5
    ),
    rep = rep(
      rep(rep(rep(1:3, each = length(mask)), 2), 5)
    ),
    method = rep(
      c(
        "Kinnex (Bambu)",
        "Kinnex (Isoquant)",
        "Kinnex (Oarfish)",
        "Kinnex (Kallisto)",
        "Kinnex (Salmon)"
      ),
      each = length(mask) * 6
    ),
    condition = rep(
      rep(c("Day 0", "Day 5"), each = length(mask) * 3),
      5
    )
  )

  set.seed(42)

  illumina_absolute_transcript <- plt_frame %>%
    filter(condition == "Day 0" & method == "Kinnex (Oarfish)") %>%
    ggplot(aes(x = illumina_abundance, y = observed_abundance, shape = as.factor(rep))) +
    geom_point(size = 4, alpha = .7) +
    geom_density_2d(color = "lightblue") +
    labs(x = "Illumina (logCPM)", y = "Kinnex (logCPM)", shape = "Replicate") +
    ggpubfigs::theme_big_simple() +
    stat_cor(aes(label = ..r.label.., colour = NULL, group = 1), show.legend = FALSE, method = "pearson", size = 6, label.x.npc = c(0.025), label.y.npc = c(1)) +
    scale_color_manual(values = ggpubfigs::friendly_pals$zesty_four[c(1, 3)]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Absolute quantification (transcript)")

  all_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("salmon_illumina", "oarfish"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/gencode/", rep, "/30000000.0/gene_counts_formatted.tsv")
        ) %>%
          arrange(desc(gene_id))
      })
    }
  ))

  prepped_counts <- lapply(
    all_counts,
    function(x) {
      x <- DGEList(x[, -1], genes = x[, 1], group = c(rep("day0", 3), rep("day5", 3)))
      cpm <- edgeR::cpm(x, log = FALSE)
      return(cpm)
    }
  )

  mask <- intersect(
    which(apply(prepped_counts[[1]], 1, function(x) sum(x > 1.0) >= 6)),
    Reduce(union, lapply(prepped_counts[-1], function(counts) {
      which(apply(counts[, 1:3], 1, function(x) (sum(x > 1.0) >= 3)) | apply(counts[, 4:6], 1, function(x) (sum(x > 1.0) >= 3)))
    }))
  )

  prepped_counts <- lapply(
    all_counts,
    function(x) {
      x <- DGEList(x[mask, -1], genes = x[mask, 1], group = c(rep("day0", 3), rep("day5", 3)))
      cpm <- edgeR::cpm(x, log = TRUE)
      return(cpm)
    }
  )



  plt_frame <- data.frame(
    log_fc = unlist(lapply(prepped_counts[-1], function(counts) {
      c(
        counts[, 1] - counts[, 4],
        counts[, 2] - counts[, 5],
        counts[, 3] - counts[, 6]
      )
    })),
    illumina_fc = rep(
      c(
        prepped_counts[[1]][, 1] - prepped_counts[[1]][, 4],
        prepped_counts[[1]][, 2] - prepped_counts[[1]][, 5],
        prepped_counts[[1]][, 3] - prepped_counts[[1]][, 6]
      ),
      1
    ),
    rep = rep(rep(1:3, each = length(mask)), 1),
    method = rep(
      c(
        "Kinnex (Oarfish)"
      ),
      each = length(mask) * 3
    )
  )


  set.seed(config$seed)

  illumina_relative_gene <- plt_frame %>%
    ggplot(aes(x = illumina_fc, y = log_fc, shape = as.factor(rep))) +
    geom_point(size = 4, alpha = .7) +
    labs(x = "Illumina (logFC)", y = "Kinnex (logFC)", shape = "Replicate") +
    ggpubfigs::theme_big_simple() +
    stat_cor(aes(label = ..r.label.., colour = NULL, group = 1), show.legend = FALSE, method = "pearson", size = 6, label.x.npc = c(0.025), label.y.npc = c(1)) +
    geom_density_2d(color = "lightblue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Relative quantification (gene)")


  plt_frame <- data.frame(
    transcript_id = all_counts[[1]]$gene_id[mask],
    observed_abundance = unlist(lapply(prepped_counts[-1], function(counts) {
      c(
        counts[, 1],
        counts[, 2],
        counts[, 3],
        counts[, 4],
        counts[, 5],
        counts[, 6]
      )
    })),
    illumina_abundance = rep(
      c(
        prepped_counts[[1]][, 1],
        prepped_counts[[1]][, 2],
        prepped_counts[[1]][, 3],
        prepped_counts[[1]][, 4],
        prepped_counts[[1]][, 5],
        prepped_counts[[1]][, 6]
      ),
      5
    ),
    rep = rep(
      rep(rep(rep(1:3, each = length(mask)), 2), 5)
    ),
    method = rep(
      c(
        "Kinnex (Bambu)",
        "Kinnex (Isoquant)",
        "Kinnex (Oarfish)",
        "Kinnex (Kallisto)",
        "Kinnex (Salmon)"
      ),
      each = length(mask) * 6
    ),
    condition = rep(
      rep(c("Day 0", "Day 5"), each = length(mask) * 3),
      5
    )
  )

  illumina_absolute_gene <- plt_frame %>%
    filter(condition == "Day 0" & method == "Kinnex (Oarfish)") %>%
    ggplot(aes(x = illumina_abundance, y = observed_abundance, shape = as.factor(rep))) +
    geom_point(size = 4, alpha = .7) +
    geom_density_2d(color = "lightblue") +
    labs(x = "Illumina (logCPM)", y = "Kinnex (logCPM)", shape = "Replicate") +
    ggpubfigs::theme_big_simple() +
    stat_cor(aes(label = ..r.label.., colour = NULL, group = 1), show.legend = FALSE, method = "pearson", size = 6, label.x.npc = c(0.025), label.y.npc = c(1)) +
    scale_color_manual(values = ggpubfigs::friendly_pals$zesty_four[c(1, 3)]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Absolute quantification (gene)")


  all_counts <- Reduce(c, lapply(
    c("oarfish", "salmon_illumina"), function(method) {
      lapply(
        c("5000000.0", "10000000.0", "20000000.0", "30000000.0"), function(depth) {
          vroom::vroom(
            paste0("results/format_quantify_subsampled/", method, "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
          )
        }
      )
    }
  ))


  transcript_gene_map <- data.frame(
    transcript_id = all_counts[[1]]$transcript_id,
    gene_id = all_counts[[1]]$gene_id
  )

  all_counts <- Reduce(c, lapply(
    c("oarfish", "run_salmon_illumina_bootstrap"), function(method) {
      lapply(
        c("5000000.0", "10000000.0", "20000000.0", "30000000.0"), function(depth) {
          if (method == "run_salmon_illumina_bootstrap") {
            s <- edgeR::catchSalmon(
              c(
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
                paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
              )
            )
            DGEList(counts = s$counts / s$annotation$Overdispersion, genes = transcript_gene_map[match(rownames(s$annotation), transcript_gene_map$transcript_id), ])
          } else {
            s <- vroom::vroom(
              paste0("results/format_quantify_subsampled/", method, "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
            )

            DGEList(
              counts = s[match(s$transcript_id, transcript_gene_map$transcript_id), -(1:2)],
              genes = transcript_gene_map[match(s$transcript_id, transcript_gene_map$transcript_id), ]
            )
          }
        }
      )
    }
  ))

  joint_isoforms <- Reduce(intersect, lapply(all_counts, function(counts) {
    which(apply(edgeR::cpm(counts$counts), 1, function(x) sum(x > 1.0) >= 3))
  }))

  all_counts_filtered <- lapply(
    all_counts, function(x) x[joint_isoforms, ]
  )

  dtu_calls <- sapply(
    all_counts,
    function(x) {
      dge <- x
      grp <- factor(c(rep(0, 3), rep(1, 3)))
      rep <- factor(rep(1:3, 2))
      mm <- model.matrix(~ grp + rep)

      keep <- filterByExpr(dge, mm)
      dge <- dge[keep, ]

      dge <- edgeR::estimateDisp(dge, mm)

      fit <- edgeR::glmQLFit(dge, mm)
      mc <- limma::makeContrasts(grp1, levels = colnames(fit$coefficients))
      ds <- diffSpliceDGE(fit, geneid = "gene_id", contrast = mc)

      tt <- topSpliceDGE(ds, number = Inf, test = "exon")

      return(tt$transcript_id[tt$FDR < 0.01])
    }
  )

  all_calls <- unique(unlist(dtu_calls))

  plt_frame <- data.frame(sapply(
    dtu_calls,
    function(x) {
      as.numeric(all_calls %in% x)
    }
  ))

  colnames(plt_frame) <- paste0(rep(c(
    "Kinnex ",
    "Illumina "
  ), each = 4), rep(c("5 M", "10 M", "20 M", "30 M"), 2))



  dtu_gencode_corrected <- upset(
    data = plt_frame,
    intersect = rev(c(
      "Kinnex 30 M",
      "Kinnex 20 M",
      "Kinnex 10 M",
      "Kinnex 5 M",
      "Illumina 30 M",
      "Illumina 20 M",
      "Illumina 10 M",
      "Illumina 5 M"
    )),
    width_ratio = 0.1, min_size = 5, wrap = TRUE, set_sizes = FALSE,
    sort_sets = FALSE,
    themes = upset_default_themes(text = element_text(size = 16)),
    base_annotations = list("Intersection size" = intersection_size(counts = FALSE)),
  ) + labs(title = "    DTU (GENCODE) (with inferential correction)") + theme(plot.title = element_text(size = 24))

  kinnex_frame <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", config$chosen_method, "/gencode/1/30000000.0/transcript_counts_formatted.tsv")
  ) %>% arrange(desc(transcript_id))

  s <- edgeR::catchSalmon(
    c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
    )
  )

  illumina_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day5-rep1` = s$counts[, 4],
    `day5-rep2` = s$counts[, 5],
    `day5-rep3` = s$counts[, 6],
    check.names = FALSE
  ) %>%
    arrange(desc(transcript_id)) %>%
    left_join(
      (
        kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
      )
    ) %>%
    dplyr::select(
      gene_id, transcript_id,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    )


  all_counts <- Reduce(c, lapply(
    c(config$chosen_method, "salmon_illumina"), function(method) {
      lapply(
        c("5000000.0", "10000000.0", "20000000.0", "30000000.0"), function(depth) {
          vroom::vroom(
            paste0("results/format_quantify_subsampled/", method, "/gencode/1/", depth, "/gene_counts_formatted.tsv")
          )
        }
      )
    }
  ))

  joint_isoforms <- Reduce(intersect, lapply(all_counts, function(counts) {
    which(apply(edgeR::cpm(counts[, -(1:2)]), 1, function(x) sum(x > 1.0) >= 3))
  }))

  all_counts_filtered <- lapply(
    all_counts, function(x) x[joint_isoforms, ]
  )

  dge_calls <- sapply(
    all_counts,
    function(x) {
      dge <- edgeR::DGEList(counts = x[, -1], genes = x[, 1], group = c(0, 0, 0, 1, 1, 1))


      grp <- sapply(colnames(dge), function(x) strsplit(x, "-")[[1]][1])
      rep <- sapply(colnames(dge), function(x) strsplit(x, "-")[[1]][2])
      mm <- model.matrix(~ grp + rep)
      keep <- edgeR::filterByExpr(dge, mm)
      dge <- dge[keep, ]
      dge <- edgeR::estimateDisp(dge, mm)

      fit <- edgeR::glmQLFit(dge, mm)
      mc <- limma::makeContrasts(grpday5, levels = colnames(fit$coefficients))

      qlf <- edgeR::glmQLFTest(fit, contrast = mc)
      tt <- edgeR::topTags(qlf, n = Inf)
      return(tt$table$gene_id[which(tt$table$FDR < 0.01)])
    }
  )

  all_calls <- unique(unlist(dge_calls))

  plt_frame <- data.frame(sapply(
    dge_calls,
    function(x) {
      as.numeric(all_calls %in% x)
    }
  ))

  colnames(plt_frame) <- paste0(rep(c(
    "Kinnex ",
    "Illumina "
  ), each = 4), rep(c("5 M", "10 M", "20 M", "30 M"), 2))



  dge_gencode <- upset(
    data = plt_frame,
    intersect = rev(c(
      "Kinnex 30 M",
      "Kinnex 20 M",
      "Kinnex 10 M",
      "Kinnex 5 M",
      "Illumina 30 M",
      "Illumina 20 M",
      "Illumina 10 M",
      "Illumina 5 M"
    )),
    width_ratio = 0.1, min_size = 100, wrap = TRUE, set_sizes = FALSE,
    sort_sets = FALSE,
    themes = upset_default_themes(text = element_text(size = 16)),
    base_annotations = list("Intersection size" = intersection_size(counts = FALSE)),
  ) + labs(title = "    DGE (GENCODE)") + theme(plot.title = element_text(size = 20))

  sixth_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        illumina_absolute_gene,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        illumina_absolute_transcript,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 2,
    labels = c("K", "", ""),
    label_size = 36
  )


  kinnex_frame <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", config$chosen_method, "/gencode/1/30000000.0/transcript_counts_formatted.tsv")
  ) %>% arrange(desc(transcript_id))

  s <- edgeR::catchSalmon(
    c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep1"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep2"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day5-rep3")
    )
  )

  illumina_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day5-rep1` = s$counts[, 4],
    `day5-rep2` = s$counts[, 5],
    `day5-rep3` = s$counts[, 6],
    check.names = FALSE
  ) %>%
    arrange(desc(transcript_id)) %>%
    left_join(
      (
        kinnex_frame %>% dplyr::select(transcript_id, gene_id) %>% distinct()
      )
    ) %>%
    dplyr::select(
      gene_id, transcript_id,
      `day0-rep1`,
      `day0-rep2`,
      `day0-rep3`,
      `day5-rep1`,
      `day5-rep2`,
      `day5-rep3`
    )


  illumina_multinomial <- illumina_frame[, 1:5] %>% arrange(transcript_id)
  kinnex_multinomial <- kinnex_frame[, 1:5] %>% arrange(transcript_id)

  ill_mask <- apply(edgeR::cpm(illumina_multinomial[, -(1:2)]), 1, function(x) all(x >= 1.0))
  kinnex_mask <- apply(edgeR::cpm(kinnex_multinomial[, -(1:2)]), 1, function(x) all(x >= 1.0))
  joint_mask <- ill_mask | kinnex_mask

  illumina_multinomial_filtered <- data.frame(illumina_multinomial[joint_mask, ], check.names = FALSE)
  kinnex_multinomial_filtered <- data.frame(kinnex_multinomial[joint_mask, ], check.names = FALSE)

  joint_genes <- intersect(illumina_multinomial_filtered$gene_id, kinnex_multinomial_filtered$gene_id)

  illumina_multinomial_filtered <- illumina_multinomial_filtered %>% filter(gene_id %in% joint_genes)
  kinnex_multinomial_filtered <- kinnex_multinomial_filtered %>% filter(gene_id %in% joint_genes)

  multi_genes <- names(which(table(illumina_multinomial_filtered$gene_id) > 1))

  illumina_cpm_rounded <- ((illumina_multinomial_filtered[, 3:5]))
  kinnex_cpm_rounded <- ((kinnex_multinomial_filtered[, 3:5]))

  dge <- DGEList(counts = cbind(as.matrix(illumina_cpm_rounded), as.matrix(kinnex_cpm_rounded)), genes = illumina_multinomial_filtered[, 1:2])

  grp <- factor(c(rep(0, 3), rep(1, 3)))
  rep <- factor(rep(1:3, 2))
  mm <- model.matrix(~ grp + rep)

  dge <- edgeR::estimateDisp(dge, mm)

  fit <- edgeR::glmQLFit(dge, mm)
  mc <- limma::makeContrasts(grp1, levels = colnames(fit$coefficients))
  ds <- diffSpliceDGE(fit, geneid = "gene_id", contrast = mc)

  tt <- topSpliceDGE(ds, number = Inf)

  tt_gene <- topSpliceDGE(ds, number = Inf, test = "gene")

  tt_exon <- topSpliceDGE(ds, number = Inf, test = "exon") %>%
    group_by(gene_id) %>%
    summarise(FDR = mean(FDR))

  p_values <- tt$FDR[match(multi_genes, tt$gene_id)]

  coldata_ill <- data.frame(
    files = c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3/quant.sf")
    ),
    names = c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3"
    )
  )

  hm <- edgeR::catchSalmon(c(
    paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep1"),
    paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep2"),
    paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/", depth, "/day0-rep3")
  ))

  se <- tximeta(coldata_ill)
  illumina_se <- fishpond::computeInfRV(se)


  kinnex_detected_original <- kinnex_multinomial$transcript_id[kinnex_mask]
  illumina_detected_original <- illumina_multinomial$transcript_id[ill_mask]


  p_value_na_mask <- !is.na(p_values)

  illumina_cpm_rounded_actual <- edgeR::cpm(illumina_cpm_rounded)

  illumina_gene_counts <- vroom::vroom(
    paste0("results/format_quantify_subsampled/", "salmon_illumina", "/gencode/1/", depth, "/gene_counts_formatted.tsv")
  )

  plt_frame <- data.frame(
    gene_id = multi_genes,
    p_value = p_values[p_value_na_mask],
    n_transcripts = sapply(multi_genes, function(gene) sum(illumina_multinomial_filtered$gene_id == gene))[p_value_na_mask],
    ill_detected_transcripts_per_gene = sapply(
      multi_genes,
      function(gene) {
        sum(illumina_multinomial_filtered$transcript_id[which(illumina_multinomial_filtered$gene_id == gene)] %in% illumina_detected_original)
      }
    )[p_value_na_mask],
    kinnex_detected_transcripts_per_gene = sapply(
      multi_genes,
      function(gene) {
        sum(kinnex_multinomial_filtered$transcript_id[which(kinnex_multinomial_filtered$gene_id == gene)] %in% kinnex_detected_original)
      }
    )[p_value_na_mask],
    kinnex_transcripts = sapply(
      multi_genes,
      function(gene) {
        (paste0(sort(kinnex_detected_original[which(kinnex_detected_original %in% kinnex_multinomial_filtered$transcript_id[which(kinnex_multinomial_filtered$gene_id == gene)])]), collapse = ","))
      }
    )[p_value_na_mask],
    illumina_transcripts = sapply(
      multi_genes,
      function(gene) {
        (paste0(sort(illumina_detected_original[which(illumina_detected_original %in% illumina_multinomial_filtered$transcript_id[which(illumina_multinomial_filtered$gene_id == gene)])]), collapse = ","))
      }
    )[p_value_na_mask],
    mean_inf_rv = sapply(
      multi_genes,
      function(gene) {
        mean(rowData(illumina_se)$meanInfRV[which(rownames(illumina_se) %in% illumina_multinomial_filtered$transcript_id[which(illumina_multinomial_filtered$gene_id == gene)])])
      }
    )[p_value_na_mask],
    overdispersion = sapply(
      multi_genes,
      function(gene) {
        mean(hm$annotation$Overdispersion[which(rownames(hm$annotation) %in% illumina_multinomial_filtered$transcript_id[which(illumina_multinomial_filtered$gene_id == gene)])])
      }
    )[p_value_na_mask],
    cpm = sapply(
      multi_genes,
      function(gene) {
        mean(as.numeric(illumina_gene_counts[which(illumina_gene_counts$gene_id == gene), 2:4]))
      }
    )[p_value_na_mask]
  )

  plt_frame$category <- ifelse(
    plt_frame$ill_detected_transcripts_per_gene > plt_frame$kinnex_detected_transcripts_per_gene,
    "Illumina overdetected",
    ifelse(plt_frame$ill_detected_transcripts_per_gene < plt_frame$kinnex_detected_transcripts_per_gene,
      "Kinnex overdetected",
      "Equal detection"
    )
  )

  plt_frame$n_transcripts <- ifelse(
    plt_frame$n_transcripts >= 5,
    "5+",
    as.character(plt_frame$n_transcripts)
  )

  split_kinnex_transcripts <- strsplit(plt_frame$kinnex_transcripts, ",")
  split_illumina_transcripts <- strsplit(plt_frame$illumina_transcripts, ",")

  hm <- ifelse(
    sapply(1:nrow(plt_frame), function(ix) all(split_kinnex_transcripts[[ix]] %in% split_illumina_transcripts[[ix]]) & all(split_illumina_transcripts[[ix]] %in% split_kinnex_transcripts[[ix]])),
    "Illumina = Kinnex",
    ifelse(
      sapply(1:nrow(plt_frame), function(ix) all(split_kinnex_transcripts[[ix]] %in% split_illumina_transcripts[[ix]])),
      "Kinnex in Illumina",
      ifelse(
        sapply(1:nrow(plt_frame), function(ix) all(split_illumina_transcripts[[ix]] %in% split_kinnex_transcripts[[ix]])),
        "Illumina in Kinnex",
        "Each subset"
      )
    )
  )

  plt_frame$better_category <- hm

  plt_frame$n_transcripts <- factor(
    plt_frame$n_transcripts,
    levels = c("2", "3", "4", "5+")
  )

  multinomial_overview <- plt_frame %>%
    filter(ill_detected_transcripts_per_gene > 0 & kinnex_detected_transcripts_per_gene > 0) %>%
    ggplot(aes(x = n_transcripts, y = -log(p_value, base = 10))) +
    geom_jitter() +
    theme_big_simple() +
    labs(x = "Transcripts per gene (GENCODE)") +
    labs(y = "-log10(q) (Mult. disc.)") +
    facet_wrap(~category) +
    geom_hline(yintercept = 5, lty = 2, color = "red", linewidth = 1)


  plt_frame$mean_factor <- ifelse(
    plt_frame$mean_inf_rv <= 5,
    "Mean inferential variability <= 5",
    "Mean inferential variability > 5"
  )

  supp_genes <- plt_frame %>%
    filter(gene_id %in% c("ENSG00000114735.10", "ENSG00000124151.19", "ENSG00000171848.16"))
  
  supp_genes$gene_name <- ifelse(
    supp_genes$gene_id == "ENSG00000114735.10",
    "HEMK1",
    ifelse(
      supp_genes$gene_id == "ENSG00000124151.19",
      "NCOA3",
      "RRM2"
    )
  )

  comparison_overview <- plt_frame %>%
    filter(ill_detected_transcripts_per_gene > 0 & kinnex_detected_transcripts_per_gene > 0) %>%
    ggplot(aes(
      y = -log(p_value, base = 10),
      #x = mean_inf_rv * (cpm + 0.01),
      x = mean_inf_rv,
      color = mean_factor
    )) +
    geom_point(size = 3) +
    theme_big_simple() +
    labs(y = "-log10(q) (Mult. disc.)", color = "") +
    #labs(x = "Mean inferential variability per gene * CPM (Illumina)") +
    labs(x = "Mean inferential variability per gene (Illumina)") +
    facet_wrap(~category) +
    scale_x_log10() +
    stat_cor(aes(label = ..r.label..), method = "pearson", size = 6, show.legend = FALSE) +
    geom_hline(yintercept = 5, lty = 2, color = "red", linewidth = 1) +
    geom_label_repel(
      data = supp_genes,
      aes(group = 1, label = gene_name, fill = NULL),
      force = 2,
      size = 7,
      nudge_y = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      "#56B4E9",
      "#CC79A7"
    ))

  first_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        ggdraw(patchwork::patchGrob(dge_gencode)),
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
                       cowplot::plot_grid(NULL,
                                          ggdraw(patchwork::patchGrob(dte_gencode)),
                                          ncol = 2,
                                          rel_widths = c(
                                            config$side_difference_large, (1 - config$side_difference_large)
                                          )
                       ),
                       nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        flips,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        flip_track,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 4,
    label_size = 36,
    labels = c("A", "B", "C", "D"),
    rel_widths = c(1, 1, 1.25, 2)
  )


  second_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
                       cowplot::plot_grid(NULL,
                                          inf_var,
                                          ncol = 2,
                                          rel_widths = c(
                                            config$side_difference_large, (1 - config$side_difference_large)
                                          )
                       ),
                       nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    #cowplot::plot_grid(NULL,
    #  cowplot::plot_grid(NULL,
    #    ggdraw(patchwork::patchGrob(dte_gencode)),
    #    ncol = 2,
    #    rel_widths = c(
    #      config$side_difference_large, (1 - config$side_difference_large)
    #    )
    #  ),
    #  nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    #),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        ggdraw(patchwork::patchGrob(dte_gencode_corrected)),
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        genes_thirty_downsampled,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        transcripts_thirty_downsampled,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        division,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 5,
    labels = c("E", "F", "G", "", "H"),
    label_size = 36,
    rel_widths = c(1, 1, 1, 1, 1.25)
  )



  third_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        (ggdraw() +
          draw_image(magick::image_read_pdf(config$division_track_adapt,
            density = config$pdf_density
          ))
        ),
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        comparison_overview,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 2,
    labels = c("I", "J"),
    label_size = 36,
    rel_widths = c(1.5, 2)
  )

  fourth_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        multinomial_overview,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        transcript_reasoning,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        overall_counts,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 3,
    labels = c("K", "L", "M"),
    label_size = 36,
    rel_widths = c(2, 1, 2)
  )

  fifth_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        illumina_relative_gene,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        illumina_relative_transcript,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        illumina_absolute_gene,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        illumina_absolute_transcript,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    rel_widths = c(1, 1, 1, 1),
    ncol = 4,
    labels = c("N", "", "O"),
    label_size = 36
  )



  final <- cowplot::plot_grid(
    first_row,
    second_row,
    third_row,
    fourth_row,
    fifth_row,
    nrow = 5,
    rel_heights = c(1, 1, 1, 1)
  )

  ggsave(output_path, final,
    dpi = 300,
    width = 26,
    height = 26
  )
}


log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(ComplexUpset)
  library(magick)
  library(ggpubfigs)
  library(dplyr)
  library(cowplot)
  library(scales)
  library(edgeR)
  library(tidyr)
  library(readr)
  library(tibble)
  library(pheatmap)
  library(tximeta)
  library(SingleCellExperiment)
  library(splines)
  library(rtracklayer)
  library(RColorBrewer)
  library(ggrepel)
  library(ggplotify)
  library(rcartocolor)
  library(ggpubr)
  library(Biostrings)
  library(tximeta)
  library(tximport)
})

sessionInfo()
status <- plot_figure_02(
  output_path = snakemake@output[[1]],
  depth = snakemake@params[["depth"]]
)

sink()
sink()

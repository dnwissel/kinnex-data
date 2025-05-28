plot_figure_S02 <- function(output_path, depth) {
  config <- rjson::fromJSON(file = "config/config.json")
  params <- rjson::fromJSON(file = "config/params.json")

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


  plt_frame <- data.frame()

  for (depth in c("5000000.0", "10000000.0", "20000000.0", "30000000.0")) {
    for (method in c("bambu", "isoquant", "salmon", "kallisto_long", "oarfish", "salmon_illumina")) {
      for (type in c("gencode")) {
        for (rep in c(
          "day0-rep1", "day0-rep2", "day0-rep3",
          "day5-rep1", "day5-rep2", "day5-rep3"
        )) {
          plt_frame <- rbind(
            plt_frame,
            data.frame(
              method = method,
              type = type,
              rep = rep,
              depth = as.numeric(depth),
              time = vroom::vroom(paste0(
                "benchmarks/quantify_subsampled/run_",
                method, "/", type, "/1/", depth, "/", rep, ".txt"
              ))$s,
              mem = vroom::vroom(paste0(
                "benchmarks/quantify_subsampled/run_",
                method, "/", type, "/1/", depth, "/", rep, ".txt"
              ))$max_uss
            )
          )
        }
      }
    }
  }

  plt_frame$method <- ifelse(
    plt_frame$method == "bambu",
    "Kinnex (Bambu)",
    ifelse(
      plt_frame$method == "isoquant",
      "Kinnex (Isoquant)",
      ifelse(
        plt_frame$method == "kallisto_long",
        "Kinnex (lr-kallisto)",
        ifelse(
          plt_frame$method == "oarfish",
          "Kinnex (Oarfish)",
          ifelse(
            plt_frame$method == "salmon_illumina",
            "Illumina",
            "Kinnex (Salmon)"
          )
        )
      )
    )
  )

  plt_frame$method <- factor(
    plt_frame$method,
    levels = c(
      "Kinnex (Bambu)", "Kinnex (Isoquant)", "Kinnex (lr-kallisto)", "Kinnex (Oarfish)", "Kinnex (Salmon)", "Illumina"
    )
  )

  time <- plt_frame %>% ggplot(aes(x = as.factor(depth), y = time / 60, color = method, group = method)) +
    geom_point(stat = "summary", fun = mean, size = 4) +
    stat_summary(fun.y = mean, geom = "line", linewidth = 2, lty = 2) +
    stat_summary(fun.data = mean_sdl, geom = "pointrange", linewidth = 1) +
    theme_big_simple() +
    labs(fill = "", y = "Runtime (m)", x = "Subsampled reads", color = "") +
    scale_color_manual(values = ggpubfigs::friendly_pals$ito_seven)

  memory <- plt_frame %>% ggplot(aes(x = as.factor(depth), y = mem / 1000, color = method, group = method)) +
    geom_point(stat = "summary", fun = mean, size = 4) +
    stat_summary(fun.y = mean, geom = "line", linewidth = 2, lty = 2) +
    stat_summary(fun.data = mean_sdl, geom = "pointrange", linewidth = 1) +
    theme_big_simple() +
    labs(fill = "", y = "Memory (USS GB)", x = "Subsampled reads", color = "") +
    scale_color_manual(values = ggpubfigs::friendly_pals$ito_seven)

  overview <- ggdraw() + draw_image(magick::image_read_pdf(config$quant_overview, density = config$pdf_density))

  exemplary_sirv_counts <- vroom::vroom(
    "results/format_quantify/salmon_illumina/sirv/transcript_counts_formatted.tsv"
  )
  short_sirvs <- grep("SIRV[0-9][0-9][0-9]$", exemplary_sirv_counts$transcript_id, value = TRUE)

  annot <- data.frame(read_excel(config$sirv_design_overview), check.names = FALSE)

  annot_df <- data.frame(
    transcript_id = annot$...6[-(1:3)][1:100],
    mw = as.numeric(annot$...8[-(1:3)][1:100]),
    e0 = as.numeric(annot$...10[-(1:3)][1:100]),
    e1 = as.numeric(annot$...11[-(1:3)][1:100]),
    e2 = as.numeric(annot$...12[-(1:3)][1:100])
  ) %>%
    drop_na() %>%
    pivot_longer(!transcript_id, names_to = "spikein", values_to = "molarity")

  annot_df_mw <- data.frame(
    transcript_id = annot$...6[-(1:3)][1:100],
    mw = as.numeric(annot$...8[-(1:3)][1:100])
  ) %>% drop_na()

  gtf_helper <- data.frame(import(config$sirv_gtf))
  sorted_transcript_length <- gtf_helper %>%
    group_by(transcript_id) %>%
    dplyr::summarise(length = sum(width)) %>%
    arrange(desc(transcript_id)) %>%
    pull(length)

  annotation_frame <- data.frame(
    transcript_id = gtf_helper %>% group_by(transcript_id) %>% dplyr::summarise(length = sum(width)) %>% arrange(desc(transcript_id)) %>% pull(transcript_id),
    len = sorted_transcript_length,
    molarity_fold_change = annot_df %>% filter(spikein == "e2") %>% arrange(desc(transcript_id)) %>% pull(molarity) / annot_df %>%
      filter(spikein == "e1") %>%
      arrange(desc(transcript_id)) %>%
      pull(molarity)
  )

  annotation_frame_relative <- data.frame(
    transcript_id = gtf_helper %>% group_by(transcript_id) %>% dplyr::summarise(length = sum(width)) %>% arrange(desc(transcript_id)) %>% pull(transcript_id),
    len = sorted_transcript_length,
    molarity_fold_change = annot_df %>% filter(spikein == "e1") %>% arrange(desc(transcript_id)) %>% pull(molarity) / annot_df %>%
      filter(spikein == "e2") %>%
      arrange(desc(transcript_id)) %>%
      pull(molarity)
  )
  annotation_frame_relative$molarity_fold_change_factor <- log(annotation_frame_relative$molarity_fold_change, base = 2)

  all_counts <- Reduce(c, lapply(
    1:5,
    function(rep) {
      lapply(c("salmon_illumina", "bambu", "isoquant", "oarfish", "kallisto", "salmon"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/sirv/", rep, "/1000000.0/transcript_counts_formatted.tsv")
        ) %>%
          filter(transcript_id %in% short_sirvs) %>%
          arrange(desc(transcript_id))
      })
    }
  ))

  au_prcs <- unlist(lapply(
    all_counts, function(counts) {
      lapply(
        c(0.01, 0.05, 0.1), function(cutoff) {
          dge <- edgeR::DGEList(counts = counts[, -(1:2)], genes = counts[, 2], group = c(0, 0, 0, 1, 1, 1))
          grp <- sapply(colnames(dge), function(x) strsplit(x, "-")[[1]][1])
          rep <- sapply(colnames(dge), function(x) strsplit(x, "-")[[1]][2])
          mm <- model.matrix(~ grp + rep)
          dge <- edgeR::estimateDisp(dge, mm)
          fit <- edgeR::glmQLFit(dge, mm)
          mc <- limma::makeContrasts(grpday5, levels = colnames(fit$coefficients))
          qlf <- edgeR::glmQLFTest(fit, contrast = mc)
          tt <- data.frame(edgeR::topTags(qlf, n = Inf)) %>% arrange(desc(transcript_id))
          truth <- annotation_frame_relative %>%
            dplyr::select(transcript_id, molarity_fold_change_factor) %>%
            distinct() %>%
            arrange(desc(transcript_id)) %>%
            mutate() %>%
            pull(molarity_fold_change_factor) %>%
            sign()
          prediction <- as.vector(decideTests(qlf, p.value = cutoff))
          positive_mask <- abs(truth) > 0
          tpr <- sum(
            truth[positive_mask] == (-1 * prediction[positive_mask])
          ) / sum(positive_mask)
          fdr <- (sum(
            abs(prediction[!positive_mask])
          ) + sum(abs(prediction[positive_mask][truth[positive_mask] == (prediction[positive_mask])]))) / sum(abs(prediction))
          list(tpr, fdr)
        }
      )
    }
  ))

  plt_frame <- data.frame(
    tpr = au_prcs[seq.int(1, length(au_prcs), 2)],
    fdr = au_prcs[seq.int(2, length(au_prcs), 2)],
    method = rep(c("Illumina", "Kinnex (Bambu)", "Kinnex (Isoquant)", "Kinnex (Oarfish)", "Kinnex (lr-kallisto)", "Kinnex (Salmon)"), each = 3),
    rep = rep(rep(1:5, each = 3), each = 6),
    cutoff = rep(rep(c(0.01, 0.05, 0.1), 6), 5)
  )

  plt_frame$method <- factor(
    plt_frame$method,
    levels = c(
      "Kinnex (Bambu)", "Kinnex (Isoquant)", "Kinnex (lr-kallisto)", "Kinnex (Oarfish)", "Kinnex (Salmon)", "Illumina"
    )
  )

  dte_sirvs <- plt_frame %>%
    group_by(cutoff, method) %>%
    summarise(tpr = mean(tpr), fdr = mean(fdr)) %>%
    ggplot(aes(x = fdr, y = tpr, color = method, group = method)) +
    geom_point(size = 4) +
    geom_path(linewidth = 1, lty = 2) +
    geom_vline(xintercept = 0.01, lty = 2) +
    geom_vline(xintercept = 0.05, lty = 2) +
    geom_vline(xintercept = 0.1, lty = 2) +
    theme_big_simple() +
    labs(fill = "", y = "TPR (SIRVs)", x = "FDR (SIRVs)", color = "") +
    scale_color_manual(values = c(ggpubfigs::friendly_pals$ito_seven[1:6])) +
    guides(colour = guide_legend(nrow = 2))

  first_row <- cowplot::plot_grid(
    cowplot::plot_grid(
      cowplot::plot_grid(NULL,
        cowplot::plot_grid(NULL,
          time + theme(legend.position = "none"),
          ncol = 2,
          rel_widths = c(
            config$side_difference_large, (1 - config$side_difference_large)
          )
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      cowplot::plot_grid(NULL,
        cowplot::plot_grid(NULL,
          memory + theme(legend.position = "none"),
          ncol = 2,
          rel_widths = c(
            config$side_difference_large, (1 - config$side_difference_large)
          )
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      cowplot::plot_grid(NULL,
        cowplot::plot_grid(NULL,
          dte_sirvs + theme(legend.position = "none"),
          ncol = 2,
          rel_widths = c(
            config$side_difference_large, (1 - config$side_difference_large)
          )
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      ncol = 3,
      rel_widths = c(0.2, 0.2, 0.25),
      labels = c("A", "", "B"),
      label_size = config$label_size
    ),
    get_plot_component(dte_sirvs, "guide-box-bottom", return_all = TRUE),
    nrow = 2,
    rel_heights = c(0.9, 0.1)
  )


  all_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("salmon_illumina", "bambu", "isoquant", "oarfish", "kallisto", "salmon"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/sirv/", rep, "/1000000.0/transcript_counts_formatted.tsv")
        ) %>%
          filter(transcript_id %in% short_sirvs) %>%
          arrange(desc(transcript_id))
      })
    }
  ))

  prepped_counts <- lapply(
    all_counts,
    function(x) {
      x <- x[x$transcript_id %in% short_sirvs, ] %>% arrange(desc(transcript_id))
      x <- DGEList(x[, -(1:2)], genes = x[, 1], group = c(rep("day0", 3), rep("day5", 3)))
      cpm <- edgeR::cpm(x, log = TRUE)
      return(cpm)
    }
  )

  annotation_frame_relative <- rbind(annotation_frame_relative, annotation_frame_relative, annotation_frame_relative)
  annotation_frame_relative$sample <- c(rep(1, 69), rep(2, 69), rep(3, 69))
  annotation_frame_relative <- rbind(
    annotation_frame_relative,
    annotation_frame_relative,
    annotation_frame_relative,
    annotation_frame_relative,
    annotation_frame_relative,
    annotation_frame_relative
  )
  annotation_frame_relative$tech <- c(
    rep("Illumina", 69 * 3),
    rep("Kinnex (Bambu)", 69 * 3),
    rep("Kinnex (Isoquant)", 69 * 3),
    rep("Kinnex (Oarfish)", 69 * 3),
    rep("Kinnex (lr-kallisto)", 69 * 3),
    rep("Kinnex (Salmon)", 69 * 3)
  )


  annotation_frame_relative$fc <- unname(unlist(lapply(prepped_counts, function(x) c(x[, 1] - x[, 4], x[, 2] - x[, 5], x[, 3] - x[, 6]))))
  annotation_frame_relative$molarity_fold_change_factor <- log(annotation_frame_relative$molarity_fold_change, base = 2)

  annotation_frame_relative$tech <- factor(
    annotation_frame_relative$tech,
    levels = c(
      "Kinnex (Bambu)",
      "Kinnex (Isoquant)",
      "Kinnex (lr-kallisto)",
      "Kinnex (Oarfish)",
      "Kinnex (Salmon)",
      "Illumina"
    )
  )

  relative <- annotation_frame_relative %>%
    ggplot(aes(x = len, y = fc, color = as.factor(molarity_fold_change_factor), shape = as.factor(sample))) +
    geom_point(size = 4, alpha = .7) +
    labs(x = "Length (SIRVs)", y = "logFC (SIRVs)", fill = "", color = "Expected logFC (SIRVs)", shape = "Replicate") +
    ggpubfigs::theme_big_simple() +
    ggplot2::scale_color_manual(values = ggpubfigs::friendly_pals$zesty_four) +
    geom_hline(yintercept = -4, col = "#F5793A", lty = 2, linewidth = 1) +
    geom_hline(yintercept = -0, col = "#A95AA1", lty = 2, linewidth = 1) +
    geom_hline(yintercept = 1, col = "#85C0F9", lty = 2, linewidth = 1) +
    geom_hline(yintercept = 6, col = "#0F2080", lty = 2, linewidth = 1) +
    geom_vline(xintercept = 1250, lty = 2, linewidth = 1, color = "red") +
    stat_cor(aes(x = as.numeric(molarity_fold_change_factor), color = NULL, y = fc, label = ..r.label.., group = 1), method = "pearson", show.legend = FALSE, label.x.npc = c(0.025), label.y.npc = c(1), size = 6) +
    facet_wrap(~tech, nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



  second_row <- cowplot::plot_grid(cowplot::plot_grid(
    cowplot::plot_grid(NULL, cowplot::plot_grid(NULL, (ggdraw() + draw_image(magick::image_read_pdf(config$fig_02_track, density = 800))), NULL, ncol = 3, rel_widths = c(0.025, 0.975, 0.025)), NULL, nrow = 3, rel_heights = c(0.05, 0.975, 0.01))
  ), nrow = 1, labels = c("C"), label_size = config$label_size)

  third_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL, relative, ncol = 2, rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    labels = c("D"),
    label_size = config$label_size
  )

  annotation_frame_absolute <- data.frame(
    transcript_id = gtf_helper %>% group_by(transcript_id) %>% dplyr::summarise(length = sum(width)) %>% arrange(desc(transcript_id)) %>% pull(transcript_id),
    len = sorted_transcript_length,
    cat = c(rep("e1", 69), rep("e2", 69)),
    molarity = c(
      annot_df %>% filter(spikein == "e1") %>% arrange(desc(transcript_id)) %>% pull(molarity),
      annot_df %>% filter(spikein == "e2") %>% arrange(desc(transcript_id)) %>% pull(molarity)
    )
  )


  annotation_frame_absolute <- rbind(
    annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute,
    annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute
  )
  annotation_frame_absolute$sample <- paste0(
    rep(rep(paste0("day", c(0, 5)), each = 69), 3),
    rep(c("-rep1", "-rep2", "-rep3"), each = 69 * 2)
  )


  annotation_frame_absolute <- rbind(
    annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute,
    annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute
  )


  all_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("salmon_illumina", "bambu", "isoquant", "oarfish", "kallisto", "salmon"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/sirv/", rep, "/1000000.0/transcript_counts_formatted.tsv")
        ) %>%
          filter(transcript_id %in% short_sirvs) %>%
          arrange(desc(transcript_id))
      })
    }
  ))

  prepped_counts <- lapply(
    all_counts,
    function(x) {
      x <- x[x$transcript_id %in% short_sirvs, ] %>% arrange(desc(transcript_id))
      x <- DGEList(x[, -(1:2)], genes = x[, 1], group = c(rep("day0", 3), rep("day5", 3)))
      cpm <- edgeR::cpm(x, log = TRUE)
      return(cpm)
    }
  )

  annotation_frame_absolute$tech <- c(
    rep("Illumina", 69 * 3 * 2),
    rep("Kinnex (Bambu)", 69 * 3 * 2), rep("Kinnex (Isoquant)", 69 * 3 * 2),
    rep("Kinnex (Oarfish)", 69 * 3 * 2), rep("Kinnex (lr-kallisto)", 69 * 3 * 2),
    rep("Kinnex (Salmon)", 69 * 3 * 2)
  )

  annotation_frame_absolute$log_cpm <- unlist(lapply(prepped_counts, function(x) {
    c(
      x[, 1],
      x[, 4],
      x[, 2],
      x[, 5],
      x[, 3],
      x[, 6]
    )
  }))

  annotation_frame_absolute$rep <- substr(sapply(strsplit(annotation_frame_absolute$sample, "\\-"), function(x) x[[2]]), 4, 4)

  annotation_frame_absolute$day <- substr(sapply(strsplit(annotation_frame_absolute$sample, "\\-"), function(x) x[[1]]), 4, 4)

  annotation_frame_absolute$tech <- factor(
    annotation_frame_absolute$tech,
    levels = c(
      "Kinnex (Bambu)",
      "Kinnex (Isoquant)",
      "Kinnex (lr-kallisto)",
      "Kinnex (Oarfish)",
      "Kinnex (Salmon)",
      "Illumina"
    )
  )

  absolute <- annotation_frame_absolute %>%
    left_join(annot_df_mw, by = "transcript_id") %>%
    filter(day == 0) %>%
    ggplot(aes(
      x = log(mw * molarity), y = log_cpm,
      shape = as.factor(rep)
    )) +
    geom_point(size = 4) +
    labs(x = "Expected abundance (SIRVs)", y = "logCPM (SIRVs)", fill = "", color = "Sample", shape = "Replicate") +
    theme_big_simple() +
    geom_smooth(aes(shape = NULL), method = "lm", se = FALSE, color = "lightblue") +
    scale_x_log10() +
    stat_cor(aes(label = ..r.label.., group = 1, color = NULL), method = "pearson", show.legend = FALSE, size = 6, label.x.npc = c(0.025), label.y.npc = c(1)) +
    facet_wrap(~tech, nrow = 1, scales = "free")

  fourth_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL, absolute, ncol = 2, rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    labels = c("E"),
    label_size = config$label_size
  )

  all_counts <- Reduce(c, lapply(
    c("salmon_illumina", "run_salmon_illumina_bootstrap", "bambu", "isoquant", "oarfish", "kallisto", "salmon"), function(method) {
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

  plt_frame <- data.frame(
    method = rep(
      c("Illumina (no correction)", "Illumina (with correction)", "Kinnex (Bambu)", "Kinnex (Isoquant)", "Kinnex (Oarfish)", "Kinnex (lr-kallisto)", "Kinnex (Salmon)"),
      each = 4
    ),
    depth = rep(
      c(5000000.0, 10000000.0, 20000000.0, 30000000.0), 7
    ),
    observed_depth = sapply(all_counts, function(x) {
      mean(apply(x$counts, 2, sum))
    })
  )

  plt_frame$method <- factor(
    plt_frame$method,
    levels = c(
      "Kinnex (Bambu)", "Kinnex (Isoquant)", "Kinnex (lr-kallisto)", "Kinnex (Oarfish)", "Kinnex (Salmon)", "Illumina (no correction)", "Illumina (with correction)"
    )
  )

  plt_frame$facet_fct <- ifelse(
    plt_frame$depth == 5e06,
    "Reads: 5 M",
    ifelse(
      plt_frame$depth == 1e07,
      "Reads: 10 M",
      ifelse(
        plt_frame$depth == 2e07,
        "Reads: 20 M",
        "Reads: 30 M"
      )
    )
  )

  plt_frame$facet_fct <- factor(
    plt_frame$facet_fct,
    levels = c("Reads: 5 M", "Reads: 10 M", "Reads: 20 M", "Reads: 30 M")
  )

  relative_quantification_efficiency <- plt_frame %>% ggplot(aes(x = method, y = observed_depth / depth, fill = method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~facet_fct, scales = "free", nrow = 1) +
    ggpubfigs::theme_big_simple() +
    labs(x = "", y = "Quantification efficiency") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(fill = "") +
    scale_fill_manual(values = ggpubfigs::friendly_pals$ito_seven)


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

  plt_frame <- data.frame(
    method = rep(
      c("Illumina (no correction)", "Illumina (with correction)", "Kinnex (Bambu)", "Kinnex (Isoquant)", "Kinnex (Oarfish)", "Kinnex (lr-kallisto)", "Kinnex (Salmon)"), 4
    ),
    depth = rep(
      c(5000000.0, 10000000.0, 20000000.0, 30000000.0),
      each = 7
    ),
    dte_calls = unlist(lapply(1:4, function(qx) {
      ix <- seq.int(qx, length.out = 7, by = 4)
      sapply(ix, function(x) {
        length(intersect(dte_calls[[x]], Reduce(union, dte_calls[ix[which(ix != x)]])))
      })
    }))
  )


  plt_frame$method <- factor(
    plt_frame$method,
    levels = c(
      "Kinnex (Bambu)", "Kinnex (Isoquant)", "Kinnex (lr-kallisto)", "Kinnex (Oarfish)", "Kinnex (Salmon)", "Illumina (no correction)", "Illumina (with correction)"
    )
  )
  
  plt_frame$facet_fct <- ifelse(
    plt_frame$depth == 5e06,
    "Reads: 5 M",
    ifelse(
      plt_frame$depth == 1e07,
      "Reads: 10 M",
      ifelse(
        plt_frame$depth == 2e07,
        "Reads: 20 M",
        "Reads: 30 M"
      )
    )
  )
  
  plt_frame$facet_fct <- factor(
    plt_frame$facet_fct,
    levels = c("Reads: 5 M", "Reads: 10 M", "Reads: 20 M", "Reads: 30 M")
  )

  significant_dte_calls <- plt_frame %>% ggplot(aes(x = method, y = dte_calls, fill = method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~facet_fct, scales = "free", nrow = 1) +
    ggpubfigs::theme_big_simple() +
    labs(x = "", y = "Significant DTE calls") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(fill = "") +
    scale_fill_manual(values = ggpubfigs::friendly_pals$ito_seven)


  fifth_row <- cowplot::plot_grid(
    cowplot::plot_grid(
      cowplot::plot_grid(NULL,
        cowplot::plot_grid(NULL, relative_quantification_efficiency + theme(legend.position = "none"), ncol = 2, rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      cowplot::plot_grid(NULL,
        cowplot::plot_grid(NULL, significant_dte_calls + theme(legend.position = "none"), ncol = 2, rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      labels = c("F", "G"),
      label_size = config$label_size
    ),
    get_plot_component(relative_quantification_efficiency, "guide-box-bottom", return_all = TRUE),
    nrow = 2,
    rel_heights = c(0.9, 0.2)
  )


  final <- cowplot::plot_grid(
    first_row,
    second_row,
    third_row,
    fourth_row,
    fifth_row,
    nrow = 5,
    rel_heights = c(1.5, 1, 1.75, 1.75, 1.75)
  )

  ggsave(output_path, final,
    dpi = 300,
    width = 26,
    height = 24
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
  library(ggplotify)
  library(rcartocolor)
  library(ggpubr)
  library(Biostrings)
  library(tximeta)
  library(tximport)
})

sessionInfo()
status <- plot_figure_S02(
  output_path = snakemake@output[[1]],
  depth = snakemake@params[["depth"]]
)

sink()
sink()

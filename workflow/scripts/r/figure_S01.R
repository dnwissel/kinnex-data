plot_figure_S01 <- function(output_path) {
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
    # Day 2-1
    adjustcolor(cols[config$exclude_color], alpha.f = 1),
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


  plt_frame <- data.frame(rbind(
    cbind(
      params$mapped_gencode_illumina,
      type = "Aligned reads", sort = "GENCODE", tech = "Illumina"
    ),
    cbind(
      params$mapped_sirv_illumina,
      type = "Aligned reads", sort = "SIRVs", tech = "Illumina"
    ),
    cbind(
      params$raw_reads_illumina,
      type = "Raw reads", sort = "All", tech = "Illumina"
    ),
    cbind(
      params$mapped_gencode_kinnex,
      type = "Aligned reads", sort = "GENCODE", tech = "Kinnex"
    ),
    cbind(
      params$mapped_sirv_kinnex,
      type = "Aligned reads", sort = "SIRVs", tech = "Kinnex"
    ),
    cbind(
      params$raw_reads_kinnex,
      type = "Raw reads", sort = "All", tech = "Kinnex"
    )
  ))

  plt_frame$rep <- c(
    "Day 0-1",
    "Day 0-2",
    "Day 0-3",
    "Day 1-1",
    "Day 2-1",
    "Day 3-1",
    "Day 3-2",
    "Day 4-1",
    "Day 5-1",
    "Day 5-2",
    "Day 5-3"
  )

  plt_frame$type <- factor(
    plt_frame$type,
    levels = c("Raw reads", "Aligned reads")
  )

  plt_frame$facet_fct <- paste0(
    plt_frame$sort, " (", plt_frame$type, ")"
  )

  plt_frame$facet_fct <- factor(
    plt_frame$facet_fct,
    levels = c(
      "All (Raw reads)",
      "GENCODE (Aligned reads)",
      "SIRVs (Aligned reads)"
    )
  )

  qc_overall <- plt_frame %>%
    filter(!facet_fct %in% c("SIRV (Quant)", "GENCODE (Quant)")) %>%
    ggplot(aes(
      x = tech, y = as.numeric(V1),
      fill = rep,
      label = paste0(
        round(as.numeric(V1) / 1e6, 2),
        " M"
      )
    )) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text(size = 4, position = position_stack(vjust = 0.5)) +
    facet_wrap(~facet_fct, scales = "free", nrow = 1) +
    ggpubfigs::theme_big_simple() +
    scale_fill_manual(
      values = replicate_cols
    ) +
    labs(x = "", y = "Reads", fill = "") +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      legend.margin = margin()
    ) +
    theme(legend.margin = margin(t = 0, unit = "cm")) +
    theme(
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0, "mm"),
      legend.spacing.y = unit(0, "mm")
    )

  quantified_gencode_kinnex <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method, "/gencode/transcript_counts_formatted.tsv")
  )

  s <- edgeR::catchSalmon(
    c(
      "results/quantify/run_salmon_illumina/gencode/day0-rep1",
      "results/quantify/run_salmon_illumina/gencode/day0-rep2",
      "results/quantify/run_salmon_illumina/gencode/day0-rep3",
      "results/quantify/run_salmon_illumina/gencode/day1-rep1",
      "results/quantify/run_salmon_illumina/gencode/day2-rep1",
      "results/quantify/run_salmon_illumina/gencode/day3-rep1",
      "results/quantify/run_salmon_illumina/gencode/day3-rep2",
      "results/quantify/run_salmon_illumina/gencode/day4-rep1",
      "results/quantify/run_salmon_illumina/gencode/day5-rep1",
      "results/quantify/run_salmon_illumina/gencode/day5-rep2",
      "results/quantify/run_salmon_illumina/gencode/day5-rep3"
    )
  )

  quantified_gencode_illumina <- data.frame(
    gene_id = NA,
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day1-rep1` = s$counts[, 4],
    `day2-rep1` = s$counts[, 5],
    `day3-rep1` = s$counts[, 6],
    `day3-rep2` = s$counts[, 7],
    `day4-rep1` = s$counts[, 8],
    `day5-rep1` = s$counts[, 9],
    `day5-rep2` = s$counts[, 10],
    `day5-rep3` = s$counts[, 11],
    check.names = FALSE
  )

  kinnex_cpm <- edgeR::cpm(quantified_gencode_kinnex[, -(1:2)])
  illumina_cpm <- edgeR::cpm(quantified_gencode_illumina[, -(1:2)])

  mask_kinnex <- apply(kinnex_cpm, 1, function(x) sum(x > config$cpm_threshold) >= config$cpm_sample_threshold)
  mask_illumina <- apply(illumina_cpm, 1, function(x) sum(x > config$cpm_threshold) >= config$cpm_sample_threshold)
  joint_mask <- mask_kinnex & mask_illumina


  cormat_illumina <- data.frame(round(cor(illumina_cpm[joint_mask, ], method = "spearman"), 3))


  cormat_kinnex <- data.frame(round(cor(kinnex_cpm[joint_mask, ], method = "spearman"), 3))

  anno_illumina <- data.frame(
    Condition = rep(c("Day 0-1", "Day 0-2", "Day 0-3", "Day 1-1", "Day 2-1", "Day 3-1", "Day 3-2", "Day 4-1", "Day 5-1", "Day 5-2", "Day 5-3"), 1)
  )

  anno_kinnex <- data.frame(
    Condition = rep(c("Day 0-1", "Day 0-2", "Day 0-3", "Day 1-1", "Day 2-1", "Day 3-1", "Day 3-2", "Day 4-1", "Day 5-1", "Day 5-2", "Day 5-3"), 1)
  )
  anno_colours <- list(
    `Condition` = c(
      `Day 0-1` = replicate_cols[1],
      `Day 0-2` = replicate_cols[2],
      `Day 0-3` = replicate_cols[3],
      `Day 1-1` = replicate_cols[4],
      `Day 2-1` = replicate_cols[5],
      `Day 3-1` = replicate_cols[6],
      `Day 3-2` = replicate_cols[7],
      `Day 4-1` = replicate_cols[8],
      `Day 5-1` = replicate_cols[9],
      `Day 5-2` = replicate_cols[10],
      `Day 5-3` = replicate_cols[11]
    )
  )
  rownames(anno_kinnex) <- rownames(cormat_kinnex)
  rownames(anno_illumina) <- rownames(cormat_illumina)

  colnames(cormat_kinnex) <- rownames(cormat_kinnex)
  colnames(cormat_illumina) <- rownames(cormat_illumina)

  heatmap_illumina <- as.ggplot(pheatmap(cormat_illumina,
    color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(10),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    show_rownames = FALSE,
    annotation_col = anno_illumina,
    annotation_row = anno_illumina,
    annotation_colors = anno_colours,
    scale = "none",
    display_numbers = TRUE,
    number_color = "black",
    main = "Illumina",
    fontsize = 16,
    fontsize_number = 9,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    height = 6,
    width = 10,
    legend = FALSE,
    annotation_legend = TRUE
  ))

  heatmap_kinnex <- as.ggplot(pheatmap(cormat_kinnex,
    color = colorRampPalette(brewer.pal(n = 7, name = "PuBuGn"))(10),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    show_rownames = FALSE,
    annotation_col = anno_kinnex,
    annotation_row = anno_kinnex,
    annotation_colors = anno_colours,
    scale = "none",
    display_numbers = TRUE,
    main = "Kinnex",
    number_color = "black",
    fontsize = 16,
    fontsize_number = 9,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    height = 6,
    width = 10,
    legend = FALSE,
    annotation_legend = TRUE
  ))

  first_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        qc_overall,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large,
          (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        heatmap_illumina,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large,
          (1 - config$side_difference_large)
        )
      ),
      nrow = 3, rel_heights = c(config$top_difference, (1 - config$top_difference), config$heatmap_difference)
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        heatmap_kinnex,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large,
          (1 - config$side_difference_large)
        )
      ),
      nrow = 3, rel_heights = c(config$top_difference, (1 - config$top_difference), config$heatmap_difference)
    ),
    ncol = 3,
    rel_widths = c(0.65, 0.25, 0.25),
    align = "v",
    labels = c("A", "B", NULL),
    label_size = config$label_size
  )

  plt_frame_qc <- data.frame()
  for (type in c("GENCODE")) {
    for (day in c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3",
      "day1-rep1",
      "day2-rep1",
      "day3-rep1",
      "day3-rep2",
      "day4-rep1",
      "day5-rep1",
      "day5-rep2",
      "day5-rep3"
    )) {
      for (tech in c("Illumina", "Kinnex")) {
        tmp_frame <- vroom::vroom(paste0(
          "results/qc/prepare_sampled_read_quality_frame_", tolower(tech), "/", day, "/", tolower(type), "/quality_plot_frame.tsv"
        ))

        if (tech == "Illumina") {
          plt_frame_qc <- rbind(
            plt_frame_qc,
            data.frame(
              tech = paste0(tech, " (", type, ")"),
              sample = day,
              tlen = tmp_frame$tlen,
              length = tmp_frame$length,
              quality = tmp_frame$quality,
              n_junctions = tmp_frame$n_junctions,
              edit_distance = tmp_frame$edit_distance
            )
          )
        } else {
          plt_frame_qc <- rbind(
            plt_frame_qc,
            data.frame(
              tech = paste0(tech, " (", type, ")"),
              sample = day,
              tlen = NA,
              length = tmp_frame$length,
              quality = tmp_frame$quality,
              n_junctions = tmp_frame$n_junctions,
              edit_distance = tmp_frame$edit_distance
            )
          )
        }
      }
    }
  }


  read_length <- plt_frame_qc %>%
    # slice_sample(n=10000) %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
        `day1-rep1` = "Day1-1",
        `day2-rep1` = "Day2-1",
        `day3-rep1` = "Day3-1",
        `day3-rep2` = "Day3-2",
        `day4-rep1` = "Day4-1",
        `day5-rep1` = "Day5-1",
        `day5-rep2` = "Day5-2",
        `day5-rep3` = "Day5-3"
      ),
      tech = recode(
        tech,
        `Illumina (GENCODE)` = "Illumina",
        `Kinnex (GENCODE)` = "Kinnex"
      )
    ) %>%
    # filter(tech == "Illumina") %>%
    ggplot(aes(
      x = sample, y = length, color = tech
    )) +
    scale_y_log10() +
    coord_flip() +
    # facet_wrap(~tech, nrow = 1) +
    ggpubfigs::theme_big_simple() +
    geom_violin(scale = "width") +
    labs(x = "", y = "Read length (bases)", color = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(
      values = c(
        "#ffb441", "#D8178C"
      )
    )

  read_bases <- plt_frame_qc %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
        `day1-rep1` = "Day1-1",
        `day2-rep1` = "Day2-1",
        `day3-rep1` = "Day3-1",
        `day3-rep2` = "Day3-2",
        `day4-rep1` = "Day4-1",
        `day5-rep1` = "Day5-1",
        `day5-rep2` = "Day5-2",
        `day5-rep3` = "Day5-3"
      ),
      tech = recode(
        tech,
        `Illumina (GENCODE)` = "Illumina",
        `Kinnex (GENCODE)` = "Kinnex"
      )
    ) %>%
    group_by(tech, sample) %>%
    summarise(total_bases = sum(length)) %>%
    ggplot(aes(
      x = sample, y = total_bases, fill = tech
    )) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    ggpubfigs::theme_big_simple() +
    labs(x = "", y = "Bases", fill = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(
      values = c(
        "#ffb441", "#D8178C"
      )
    )


  junctions <- plt_frame_qc %>%
    filter(tech %in% c("Illumina (GENCODE)", "Kinnex (GENCODE)")) %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
        `day1-rep1` = "Day1-1",
        `day2-rep1` = "Day2-1",
        `day3-rep1` = "Day3-1",
        `day3-rep2` = "Day3-2",
        `day4-rep1` = "Day4-1",
        `day5-rep1` = "Day5-1",
        `day5-rep2` = "Day5-2",
        `day5-rep3` = "Day5-3"
      ),
      tech = recode(
        tech,
        `Illumina (GENCODE)` = "Illumina",
        `Kinnex (GENCODE)` = "Kinnex"
      )
    ) %>%
    ggplot(aes(
      x = n_junctions, y = sample, color = tech
    )) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    ggpubfigs::theme_big_simple() +
    geom_violin(scale = "width") +
    labs(y = "", x = "Number of covered junctions", color = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(
      values = c(
        "#ffb441", "#D8178C"
      )
    )

  tlengths <- plt_frame_qc %>%
    filter(tech %in% c("Illumina (GENCODE)")) %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
        `day1-rep1` = "Day1-1",
        `day2-rep1` = "Day2-1",
        `day3-rep1` = "Day3-1",
        `day3-rep2` = "Day3-2",
        `day4-rep1` = "Day4-1",
        `day5-rep1` = "Day5-1",
        `day5-rep2` = "Day5-2",
        `day5-rep3` = "Day5-3"
      ),
      tech = recode(
        tech,
        `Illumina (GENCODE)` = "Illumina",
        `Kinnex (GENCODE)` = "Kinnex"
      )
    ) %>%
    ggplot(aes(
      x = tlen, y = sample, color = tech
    )) +
    scale_x_log10() +
    ggpubfigs::theme_big_simple() +
    geom_violin(scale = "width") +
    labs(y = "", x = "Fragment size per read", color = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(
      values = c(
        "#ffb441", "#D8178C"
      )
    )

  three_prime_bias_illumina_gencode <- data.frame(
    `day0-rep1` = params$illumina_day0_rep1,
    `day0-rep2` = params$illumina_day0_rep2,
    `day0-rep3` = params$illumina_day0_rep3,
    `day1-rep1` = params$illumina_day1_rep1,
    `day2-rep1` = params$illumina_day2_rep1,
    `day3-rep1` = params$illumina_day3_rep1,
    `day3-rep2` = params$illumina_day3_rep2,
    `day4-rep1` = params$illumina_day4_rep1,
    `day5-rep1` = params$illumina_day5_rep1,
    `day5-rep2` = params$illumina_day5_rep2,
    `day5-rep3` = params$illumina_day5_rep3
  )

  three_prime_bias_kinnex_gencode <- data.frame(
    `day0-rep1` = params$kinnex_day0_rep1,
    `day0-rep2` = params$kinnex_day0_rep2,
    `day0-rep3` = params$kinnex_day0_rep3,
    `day1-rep1` = params$kinnex_day1_rep1,
    `day2-rep1` = params$kinnex_day2_rep1,
    `day3-rep1` = params$kinnex_day3_rep1,
    `day3-rep2` = params$kinnex_day3_rep2,
    `day4-rep1` = params$kinnex_day4_rep1,
    `day5-rep1` = params$kinnex_day5_rep1,
    `day5-rep2` = params$kinnex_day5_rep2,
    `day5-rep3` = params$kinnex_day5_rep3
  )


  rbind(
    cbind(three_prime_bias_illumina_gencode, type = "Illumina", tech = "GE"),
    cbind(three_prime_bias_kinnex_gencode, type = "Kinnex", tech = "GE")
  ) -> plt_frame

  plt_frame$x <- 1:100

  plt_frame <- plt_frame %>%
    pivot_longer(
      cols = starts_with("day"),
      names_to = "day",
      values_to = "val"
    ) %>%
    group_by(day, tech, type) %>%
    arrange(x, .by_group = TRUE)

  plt_frame$day_formatted <- paste0("Day ", sapply(
    strsplit(plt_frame$day, "\\."),
    function(x) substr(x[[1]], nchar(x[[1]]), nchar(x[[1]]))
  ))

  plt_frame$rep_formatted <- sapply(
    strsplit(plt_frame$day, "\\."),
    function(x) substr(x[[2]], nchar(x[[2]]), nchar(x[[2]]))
  )

  plt_frame$day_formatted <- paste0(plt_frame$day_formatted, "-", plt_frame$rep_formatted)
  three_prime_bias <- plt_frame %>%
    filter(tech == "GE") %>%
    ggplot(aes(x = as.numeric(x), y = as.numeric(val), col = day_formatted)) +
    geom_path(linewidth = 2) +
    facet_wrap(~type) +
    labs(x = "Gene body percentile (5' to 3')", y = "Norm. coverage") +
    ggpubfigs::theme_big_simple() +
    scale_color_manual(
      values = replicate_cols
    ) +
    labs(color = "") +
    guides(col = guide_legend(nrow = 3, byrow = TRUE))


  read_qual <- plt_frame_qc %>%
    filter(tech %in% c("Illumina (GENCODE)", "Kinnex (GENCODE)")) %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
        `day1-rep1` = "Day1-1",
        `day2-rep1` = "Day2-1",
        `day3-rep1` = "Day3-1",
        `day3-rep2` = "Day3-2",
        `day4-rep1` = "Day4-1",
        `day5-rep1` = "Day5-1",
        `day5-rep2` = "Day5-2",
        `day5-rep3` = "Day5-3"
      ),
      tech = recode(
        tech,
        `Illumina (GENCODE)` = "Illumina",
        `Kinnex (GENCODE)` = "Kinnex"
      )
    ) %>%
    ggplot(aes(
      x = sample, y = quality, color = tech
    )) +
    geom_violin(scale = "width") +
    scale_y_log10() +
    coord_flip() +
    ggpubfigs::theme_big_simple() +
    labs(x = "", y = "Read quality (average phred)", color = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(
      values = c(
        "#ffb441", "#D8178C"
      )
    )

  edit_distance <- plt_frame_qc %>%
    filter(tech %in% c("Illumina (GENCODE)", "Kinnex (GENCODE)")) %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
        `day1-rep1` = "Day1-1",
        `day2-rep1` = "Day2-1",
        `day3-rep1` = "Day3-1",
        `day3-rep2` = "Day3-2",
        `day4-rep1` = "Day4-1",
        `day5-rep1` = "Day5-1",
        `day5-rep2` = "Day5-2",
        `day5-rep3` = "Day5-3"
      ),
      tech = recode(
        tech,
        `Illumina (GENCODE)` = "Illumina",
        `Kinnex (GENCODE)` = "Kinnex"
      )
    ) %>%
    ggplot(aes(
      x = as.numeric(edit_distance), y = sample, color = tech
    )) +
    geom_violin(scale = "width") +
    coord_cartesian(xlim = c(0, 0.05)) +
    ggpubfigs::theme_big_simple() +
    labs(y = "", x = "Relative edit distance", color = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(
      values = c(
        "#ffb441", "#D8178C"
      )
    )

  second_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL, read_length, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL, read_bases, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL, junctions, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL, tlengths, ncol = 2, rel_widths = c(0.05, 0.975)),
    nrow = 1, ncol = 4,
    align = "hv",
    label_size = 36,
    labels = c("C", "D", "E", "F")
  )


  quantified_gencode_kinnex <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method, "/gencode/transcript_counts_formatted.tsv")
  )

  s <- edgeR::catchSalmon(
    c(
      "results/quantify/run_salmon_illumina/gencode/day0-rep1",
      "results/quantify/run_salmon_illumina/gencode/day0-rep2",
      "results/quantify/run_salmon_illumina/gencode/day0-rep3",
      "results/quantify/run_salmon_illumina/gencode/day1-rep1",
      "results/quantify/run_salmon_illumina/gencode/day2-rep1",
      "results/quantify/run_salmon_illumina/gencode/day3-rep1",
      "results/quantify/run_salmon_illumina/gencode/day3-rep2",
      "results/quantify/run_salmon_illumina/gencode/day4-rep1",
      "results/quantify/run_salmon_illumina/gencode/day5-rep1",
      "results/quantify/run_salmon_illumina/gencode/day5-rep2",
      "results/quantify/run_salmon_illumina/gencode/day5-rep3"
    )
  )

  quantified_gencode_illumina <- data.frame(
    gene_id = NA,
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day1-rep1` = s$counts[, 4],
    `day2-rep1` = s$counts[, 5],
    `day3-rep1` = s$counts[, 6],
    `day3-rep2` = s$counts[, 7],
    `day4-rep1` = s$counts[, 8],
    `day5-rep1` = s$counts[, 9],
    `day5-rep2` = s$counts[, 10],
    `day5-rep3` = s$counts[, 11],
    check.names = FALSE
  )


  kinnex_dgelist <- DGEList(counts = quantified_gencode_kinnex[, -(1:2)], group = c(
    0, 0, 0,
    1,
    2,
    3, 3,
    4,
    5, 5, 5
  ))

  grp <- sapply(strsplit(colnames(kinnex_dgelist), "-", fixed = TRUE), .subset, 1)
  reps <- sapply(strsplit(colnames(kinnex_dgelist), "-", fixed = TRUE), .subset, 2)

  mm <- model.matrix(~ 0 + grp + reps)
  kinnex_dgelist <- calcNormFactors(kinnex_dgelist)
  kinnex_dgelist <- estimateDisp(kinnex_dgelist, mm)
  plotBCV(kinnex_dgelist)

  tmp <- plotMDS(kinnex_dgelist, cex = 2)
  df <- data.frame(
    Sample = rownames(tmp$distance.matrix.squared),
    MDS1 = tmp$x, MDS2 = tmp$y
  )

  ill_dgelist <- DGEList(
    counts = quantified_gencode_illumina[, -(1:2)],
    group = c(
      0, 0, 0,
      1,
      2,
      3, 3,
      4,
      5, 5, 5
    )
  )

  grp <- sapply(strsplit(colnames(ill_dgelist), "-", fixed = TRUE), .subset, 1)
  mm <- model.matrix(~ 0 + grp + reps)
  ill_dgelist <- calcNormFactors(ill_dgelist)
  ill_dgelist <- estimateDisp(ill_dgelist, mm)
  tmp <- plotMDS(ill_dgelist, cex = 2)

  design <- cbind(Intercept = 1, xaxis = tmp$x, yaxis = tmp$y)
  fit <- lmFit(ill_dgelist$counts, design)
  fit <- eBayes(fit)
  df_illumina <- data.frame(
    Sample = rownames(tmp$distance.matrix.squared),
    MDS1 = tmp$x, MDS2 = tmp$y
  )

  df_overall <- rbind(
    cbind(df, type = "Kinnex"),
    cbind(df_illumina, type = "Illumina")
  )

  df_overall$Sample <- paste0(
    "Day ", substr(sapply(strsplit(
      df_overall$Sample, "-"
    ), function(x) x[[1]]), 4, 4), "-",
    substr(sapply(strsplit(df_overall$Sample, "-"), function(x) x[[2]]), 4, 4)
  )


  mds_gencode <- cowplot::plot_grid(
    df_overall %>%
      filter(type == "Illumina") %>%
      mutate(type = "Illumina") %>%
      ggplot(aes(x = MDS1, y = MDS2, label = Sample)) +
      geom_point(size = 4) +
      geom_label_repel(min.segment.length = 0, size = 5, max.overlaps = 15) +
      labs(x = "MDS1 (21%)", y = "MDS2 (11%)") +
      facet_wrap(~type) +
      theme_big_simple(),
    df_overall %>%
      filter(type == "Kinnex") %>%
      mutate(type = "Kinnex") %>%
      ggplot(aes(x = MDS1, y = MDS2, label = Sample)) +
      geom_point(size = 4) +
      geom_label_repel(min.segment.length = 0, size = 5, max.overlaps = 15) +
      labs(x = "MDS1 (35%)", y = "MDS2 (22%)") +
      facet_wrap(~type) +
      theme_big_simple()
  )

  third_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL, read_qual, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL, edit_distance, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL, three_prime_bias, ncol = 2, rel_widths = c(0.05, 0.975)),
    nrow = 1, ncol = 3,
    rel_widths = c(1, 1, 2),
    align = "hv",
    label_size = 36,
    labels = c("G", "H", "I")
  )
  fourth_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL, mds_gencode, ncol = 2, rel_widths = c(0.05, 0.975)),
    nrow = 1, ncol = 1,
    rel_widths = c(1),
    align = "hv",
    label_size = 36,
    labels = c("J")
  )

  final <- cowplot::plot_grid(
    first_row,
    second_row,
    third_row,
    fourth_row,
    nrow = 4,
    rel_heights = c(1, 0.75, 0.75, 0.5)
  )

  ggsave(output_path, final,
    dpi = 300,
    width = 26,
    height = 22
  )
}


log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(readr)
  library(ggrepel)
  library(edgeR)
  library(pheatmap)
  library(reshape2)
  library(splines)
  library(SummarizedExperiment)
  library(rcartocolor)
  library(cowplot)
  library(rtracklayer)
  library(tidyr)
  library(ggpmisc)
  library(biomaRt)
  library(magick)
  library(scales)
  library(ggpubfigs)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplotify)
})

sessionInfo()
status <- plot_figure_S01(
  output_path = snakemake@output[[1]]
)

sink()
sink()

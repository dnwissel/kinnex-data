plot_figure_01 <- function(output_path) {
  config <- rjson::fromJSON(file = "config/config.json")
  params <- rjson::fromJSON(file = "config/params.json")

  cols <- carto_pal(config$n_colors, "Safe")
  replicate_cols <- cols[-config$exclude_color]
  replicate_cols <- c(replicate_cols[1], replicate_cols[5], replicate_cols[3], replicate_cols[4], replicate_cols[2])

  replicate_cols <- c(
    adjustcolor(replicate_cols[1], alpha.f = 1),
    adjustcolor(replicate_cols[1], alpha.f = 0.75),
    adjustcolor(replicate_cols[1], alpha.f = 0.5),
    adjustcolor(replicate_cols[2], alpha.f = 1),
    adjustcolor(replicate_cols[3], alpha.f = 1),
    adjustcolor(replicate_cols[3], alpha.f = 0.75),
    adjustcolor(replicate_cols[4], alpha.f = 1),
    adjustcolor(replicate_cols[5], alpha.f = 1),
    adjustcolor(replicate_cols[5], alpha.f = 0.75),
    adjustcolor(replicate_cols[5], alpha.f = 0.5)
  )

  first_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL, ggdraw() +
        draw_image(magick::image_read_pdf(config$exp_design_path,
          density = config$pdf_density
        )),
      ncol = 2, rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
      ),
      nrow = 2, rel_heights = c(
        config$top_difference,
        (1 - config$top_difference)
      )
    ),
    labels = c("A"), label_size = config$label_size
  )

  plt_frame <- data.frame(rbind(
    cbind(
      params$mapped_gencode_illumina[-5],
      type = "Aligned reads", sort = "GENCODE", tech = "Illumina"
    ),
    cbind(
      params$mapped_sirv_illumina[-5],
      type = "Aligned reads", sort = "SIRVs", tech = "Illumina"
    ),
    cbind(
      params$raw_reads_illumina[-5],
      type = "Raw reads", sort = "All", tech = "Illumina"
    ),
    cbind(
      params$mapped_gencode_kinnex[-5],
      type = "Aligned reads", sort = "GENCODE", tech = "Kinnex"
    ),
    cbind(
      params$mapped_sirv_kinnex[-5],
      type = "Aligned reads", sort = "SIRVs", tech = "Kinnex"
    ),
    cbind(
      params$raw_reads_kinnex[-5],
      type = "Raw reads", sort = "All", tech = "Kinnex"
    )
  ))

  plt_frame$rep <- c(
    "Day 0-1",
    "Day 0-2",
    "Day 0-3",
    "Day 1-1",
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
    ) +
    theme(axis.title.x = element_blank())

  quantified_gencode_kinnex <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method, "/gencode/transcript_counts_formatted.tsv")
  )[, -7] %>% arrange(desc(transcript_id))

  s <- edgeR::catchSalmon(
    c(
      "results/quantify/run_salmon_illumina/gencode/day0-rep1",
      "results/quantify/run_salmon_illumina/gencode/day0-rep2",
      "results/quantify/run_salmon_illumina/gencode/day0-rep3",
      "results/quantify/run_salmon_illumina/gencode/day1-rep1",
      "results/quantify/run_salmon_illumina/gencode/day3-rep1",
      "results/quantify/run_salmon_illumina/gencode/day3-rep2",
      "results/quantify/run_salmon_illumina/gencode/day4-rep1",
      "results/quantify/run_salmon_illumina/gencode/day5-rep1",
      "results/quantify/run_salmon_illumina/gencode/day5-rep2",
      "results/quantify/run_salmon_illumina/gencode/day5-rep3"
    )
  )

  quantified_gencode_illumina <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day1-rep1` = s$counts[, 4],
    `day3-rep1` = s$counts[, 5],
    `day3-rep2` = s$counts[, 6],
    `day4-rep1` = s$counts[, 7],
    `day5-rep1` = s$counts[, 8],
    `day5-rep2` = s$counts[, 9],
    `day5-rep3` = s$counts[, 10],
    check.names = FALSE
  ) %>%
    left_join(
      data.frame(
        gene_id = quantified_gencode_kinnex$gene_id,
        transcript_id = quantified_gencode_kinnex$transcript_id
      ) %>% distinct()
    ) %>%
    dplyr::select(
      gene_id,
      transcript_id,
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
    arrange(desc(transcript_id))

  kinnex_cpm <- edgeR::cpm(quantified_gencode_kinnex[, -(1:2)])
  illumina_cpm <- edgeR::cpm(quantified_gencode_illumina[, -(1:2)])

  mask_kinnex <- apply(kinnex_cpm, 1, function(x) sum(x > config$cpm_threshold) >= config$cpm_sample_threshold)
  mask_illumina <- apply(illumina_cpm, 1, function(x) sum(x > config$cpm_threshold) >= config$cpm_sample_threshold)
  joint_mask <- mask_kinnex & mask_illumina


  cormat_illumina <- data.frame(round(cor(illumina_cpm[joint_mask, ], method = "spearman"), 3))


  cormat_kinnex <- data.frame(round(cor(kinnex_cpm[joint_mask, ], method = "spearman"), 3))

  anno_illumina <- data.frame(
    Condition = rep(c("Day 0-1", "Day 0-2", "Day 0-3", "Day 1-1", "Day 3-1", "Day 3-2", "Day 4-1", "Day 5-1", "Day 5-2", "Day 5-3"), 1)
  )

  anno_kinnex <- data.frame(
    Condition = rep(c("Day 0-1", "Day 0-2", "Day 0-3", "Day 1-1", "Day 3-1", "Day 3-2", "Day 4-1", "Day 5-1", "Day 5-2", "Day 5-3"), 1)
  )
  anno_colours <- list(
    `Condition` = c(
      `Day 0-1` = replicate_cols[1],
      `Day 0-2` = replicate_cols[2],
      `Day 0-3` = replicate_cols[3],
      `Day 1-1` = replicate_cols[4],
      `Day 3-1` = replicate_cols[5],
      `Day 3-2` = replicate_cols[6],
      `Day 4-1` = replicate_cols[7],
      `Day 5-1` = replicate_cols[8],
      `Day 5-2` = replicate_cols[9],
      `Day 5-3` = replicate_cols[10]
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
    number_color = "black",
    fontsize = 16,
    fontsize_number = 9,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    height = 6,
    main = "Kinnex",
    width = 10,
    legend = FALSE,
    annotation_legend = TRUE
  ))

  second_row <- cowplot::plot_grid(
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
    labels = c("B", "J", NULL),
    label_size = config$label_size
  )


  plt_frame_qc <- data.frame()
  for (type in c("GENCODE")) {
    for (day in c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3",
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
    filter(tech %in% c("Illumina (GENCODE)", "Kinnex (GENCODE)")) %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
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
      x = sample, y = length, color = tech
    )) +
    scale_y_log10() +
    coord_flip() +
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
    filter(tech %in% c("Illumina (GENCODE)", "Kinnex (GENCODE)")) %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
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

  kinnex_cpm <- edgeR::cpm(quantified_gencode_kinnex[, -(1:2)], log = TRUE)
  illumina_cpm <- edgeR::cpm(quantified_gencode_illumina[, -(1:2)], log = TRUE)

  quantified_gencode_kinnex_gene <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method, "/gencode/gene_counts_formatted.tsv")
  )[, -6]

  quantified_gencode_illumina_gene <- vroom::vroom(
    "results/format_quantify/salmon_illumina/gencode/gene_counts_formatted.tsv"
  )[, -6]

  kinnex_cpm_gene <- edgeR::cpm(quantified_gencode_kinnex_gene[, -1])
  illumina_cpm_gene <- edgeR::cpm(quantified_gencode_illumina_gene[, -1])

  mask_kinnex_gene <- apply(kinnex_cpm_gene, 1, function(x) sum(x > config$cpm_threshold) >= config$cpm_sample_threshold)
  mask_illumina_gene <- apply(illumina_cpm_gene, 1, function(x) sum(x > config$cpm_threshold) >= config$cpm_sample_threshold)
  joint_mask_gene <- mask_kinnex_gene & mask_illumina_gene

  kinnex_cpm_gene <- edgeR::cpm(quantified_gencode_kinnex_gene[, -1], log = TRUE)
  illumina_cpm_gene <- edgeR::cpm(quantified_gencode_illumina_gene[, -1], log = TRUE)

  plt_frame <- data.frame(
    `Day0-1 (logCPM)` = c(
      illumina_cpm[joint_mask, 1],
      kinnex_cpm[joint_mask, 1],
      illumina_cpm_gene[joint_mask_gene, 1],
      kinnex_cpm_gene[joint_mask_gene, 1]
    ),
    `Day0-2 (logCPM)` = c(
      illumina_cpm[joint_mask, 2],
      kinnex_cpm[joint_mask, 2],
      illumina_cpm_gene[joint_mask_gene, 2],
      kinnex_cpm_gene[joint_mask_gene, 2]
    ),
    type = c(
      rep("Illumina (transcript)", sum(joint_mask)),
      rep("Kinnex (transcript)", sum(joint_mask)),
      rep("Illumina (gene)", sum(joint_mask_gene)),
      rep("Kinnex (gene)", sum(joint_mask_gene))
    ),
    check.names = FALSE
  )

  plt_frame$type <- factor(
    plt_frame$type,
    levels = c("Illumina (gene)", "Kinnex (gene)", "Illumina (transcript)", "Kinnex (transcript)")
  )

  replicability <- plt_frame %>% ggplot(aes(x = `Day0-1 (logCPM)`, y = `Day0-2 (logCPM)`)) +
    geom_point(alpha = 0.7) +
    stat_cor(aes(label = ..r.label..), method = "spearman", show.legend = FALSE, label.y = 15, label.x = -5, size = 6) +
    facet_wrap(~type) +
    geom_density_2d(color = "lightblue") +
    theme_big_simple()

  three_prime_bias_illumina_gencode <- data.frame(
    `day0-rep1` = params$illumina_day0_rep1,
    `day0-rep2` = params$illumina_day0_rep2,
    `day0-rep3` = params$illumina_day0_rep3,
    `day5-rep1` = params$illumina_day5_rep1,
    `day5-rep2` = params$illumina_day5_rep2,
    `day5-rep3` = params$illumina_day5_rep3
  )

  three_prime_bias_kinnex_gencode <- data.frame(
    `day0-rep1` = params$kinnex_day0_rep1,
    `day0-rep2` = params$kinnex_day0_rep2,
    `day0-rep3` = params$kinnex_day0_rep3,
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
    filter(day_formatted != "Day 2-1") %>%
    ggplot(aes(x = as.numeric(x), y = as.numeric(val), col = day_formatted)) +
    geom_path(linewidth = 2) +
    facet_wrap(~type) +
    labs(x = "Gene body percentile (5' to 3')", y = "Norm. coverage") +
    ggpubfigs::theme_big_simple() +
    scale_color_manual(
      values = replicate_cols[c(1:3, 8:10)]
    ) +
    labs(color = "") +
    guides(col = guide_legend(nrow = 2, byrow = TRUE))

  intermediate_restructure_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL, read_length, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL, read_bases, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL, junctions, ncol = 2, rel_widths = c(0.05, 0.975)),
    nrow = 1, ncol = 3,
    rel_widths = c(1, 1, 1.4),
    align = "hv",
    label_size = config$label_size,
    labels = c("C", "D", "E")
  )


  read_qual <- plt_frame_qc %>%
    filter(tech %in% c("Illumina (GENCODE)", "Kinnex (GENCODE)")) %>%
    mutate(
      sample = recode(
        sample,
        `day0-rep1` = "Day0-1",
        `day0-rep2` = "Day0-2",
        `day0-rep3` = "Day0-3",
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

  intermediate_restructure_row_two <- cowplot::plot_grid(
    cowplot::plot_grid(NULL, read_qual, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL, edit_distance, ncol = 2, rel_widths = c(0.05, 0.975)),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        three_prime_bias,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large * 2, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    nrow = 1, ncol = 3,
    rel_widths = c(1, 1, 1.4),
    label_size = config$label_size,
    labels = c("F", "G", "H")
  )


  third_row <- cowplot::plot_grid(
    cowplot::plot_grid(
      intermediate_restructure_row,
      intermediate_restructure_row_two,
      nrow = 2,
      rel_heights = c(0.5, 0.5),
      align = "v"
    ),
    NULL,
    cowplot::plot_grid(NULL, cowplot::plot_grid(NULL,
      replicability,
      ncol = 2,
      rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
    ),
    nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 3,
    rel_widths = c(0.75, 0.025, 0.4),
    labels = c("", "K"),
    label_size = config$label_size
  )

  novel_counts <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method, "/gencode/gene_counts_formatted.tsv")
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
    ), multiple = "first") %>%
    pull(gene_name)


  y <- DGEList(
    counts = novel_counts[, c(2:4, 10:12)],
    genes = data.frame(gene_id = novel_counts$gene_id, gene_name = gene_name)
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
  tt$table$gene_name[is.na(tt$table$gene_name)] <- ""
  plt_frame <- data.frame(
    coef = tt$table$logFC,
    p_val = tt$table$FDR,
    gene_id = tt$table$gene_id,
    gene_name = tt$table$gene_name
  )

  plt_frame$gene_type <- ifelse(
    plt_frame$p_val < 0.01 & plt_frame$coef > 2,
    "Up",
    ifelse(
      plt_frame$p_val < 0.01 & plt_frame$coef < -2,
      "Down",
      "NS"
    )
  )

  cols <- c("Up" = config$up_down_cols[2], "Down" = config$up_down_cols[1], "NS" = "grey")
  sizes <- c("Up" = 2, "Down" = 2, "NS" = 1)
  alphas <- c("Up" = 1, "Down" = 1, "NS" = 0.5)

  sig_il_genes <- plt_frame %>%
    filter(gene_name %in% c("NANOG", "ETV2", "KDR", "CD34", "SOX17", "SOX7"))

  # Adapted from: https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/
  vol_plot <- plt_frame %>%
    ggplot(aes(
      x = coef,
      y = -log10(p_val),
    )) +
    geom_point(aes(colour = gene_type),
      alpha = 0.7,
      shape = 16,
      size = 4
    ) +
    geom_hline(
      yintercept = -log10(0.01),
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c(log2(0.25), log2(4)),
      linetype = "dashed"
    ) +
    xlim(-17.5, 17.5) +
    geom_label_repel(
      data = sig_il_genes,
      aes(label = gene_name, fill = NULL),
      force = 2,
      size = 7,
      nudge_y = 1
    ) +
    scale_color_manual(values = cols) +
    ggpubfigs::theme_big_simple() +
    labs(x = "logFC (Day 5 vs Day 0)", y = "-log10(p)", color = "")


  novel_counts <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method, "/gencode/transcript_counts_formatted.tsv")
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
  goi <- "ENSG00000086200.17"
  w <- y$genes$gene_id == goi
  gn <- y$genes$gene_name[y$genes$gene_id == goi][1]
  fitted <- y$fitted.values
  rownames(cps) <- rownames(fitted) <- y$genes$tx_id
  colnames(cps) <- rownames(y$samples)
  md$sample_id <- paste0("day", md$day, "-", md$rep)

  m <- melt(fitted[w, ]) %>%
    setNames(c("txid", "sample_id", "fitted")) %>%
    left_join(md)
  z <- mm[, 4:5] %*% t(y$coefficients[, 4:5])

  a <- melt(z[, w]) %>% setNames(c("md_row", "tx_id", "rel_expr"))
  a$tx_id <- ifelse(
    a$tx_id == 20585,
    sapply(strsplit(y$genes$tx_id[which(w)][2], "\\."), function(x) x[[1]]),
    sapply(strsplit(y$genes$tx_id[which(w)][1], "\\."), function(x) x[[1]])
  )

  a$tx_id <- factor(
    a$tx_id,
    levels = c(
      sapply(strsplit(y$genes$tx_id[which(w)][2], "\\."), function(x) x[[1]]),
      sapply(strsplit(y$genes$tx_id[which(w)][1], "\\."), function(x) x[[1]])
    ),
    ordered = TRUE
  )

  p_final <- ggplot(a, aes(
    x = md$day[md_row], y = rel_expr,
    colour = as.factor(tx_id), group = tx_id
  )) +
    geom_point(size = 4) +
    geom_line(size = 2) +
    ggtitle(paste("Model fits: ", gn)) +
    xlab("day") +
    labs(x = "Day", y = "Rel. expression", title = "", color = "") +
    ggpubfigs::theme_big_simple() +
    scale_color_manual(values = (config$up_down_cols)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

  b <- melt(cps[w, ]) %>%
    setNames(c("tx_id", "sample_id", "cpm")) %>%
    left_join(md)
  b$tx_id <- factor(
    b$tx_id,
    levels = rev(c(y$genes$tx_id[which(w)][2], y$genes$tx_id[which(w)][1])),
    ordered = TRUE
  )

  q_final <- ggplot(b, aes(
    x = day, y = cpm + .5,
    colour = tx_id,
    shape = rep
  )) +
    geom_point(size = 4) +
    scale_y_log10() +
    ggtitle(paste("CPMs: ", gn)) +
    labs(x = "Day", y = "CPM", title = "") +
    ggpubfigs::theme_big_simple() +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_color_manual(values = (config$up_down_cols)) +
    labs(color = "", shape = "") +
    guides(colour = "none")

  spline_plot <- cowplot::plot_grid(p_final, q_final, ncol = 2)

  kinnex_dgelist <- DGEList(counts = quantified_gencode_kinnex[, -(1:2)], group = c(
    0, 0, 0,
    1,
    3, 3,
    4,
    5, 5, 5
  ))

  grp <- sapply(strsplit(colnames(kinnex_dgelist), "-", fixed = TRUE), .subset, 1)

  mm <- model.matrix(~ 0 + grp)
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
      3, 3,
      4,
      5, 5, 5
    )
  )

  grp <- sapply(strsplit(colnames(ill_dgelist), "-", fixed = TRUE), .subset, 1)
  mm <- model.matrix(~ 0 + grp)
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
      labs(x = "MDS1 (23%)", y = "MDS2 (12%)") +
      facet_wrap(~type) +
      theme_big_simple(),
    df_overall %>%
      filter(type == "Kinnex") %>%
      mutate(type = "Kinnex") %>%
      ggplot(aes(x = MDS1, y = MDS2, label = Sample)) +
      geom_point(size = 4) +
      geom_label_repel(min.segment.length = 0, size = 5, max.overlaps = 15) +
      labs(x = "MDS1 (45%)", y = "MDS2 (13%)") +
      facet_wrap(~type) +
      theme_big_simple()
  )


  novel_counts <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method_novel, "_novel/gencode/transcript_counts_formatted.tsv")
  )[, -7]


  transcriptome_coding <- import(
    "results/annotate/run_orfanage/transcriptome.gtf"
  )


  cps <- edgeR::cpm(novel_counts[, -(1:2)])[grep("Bambu", novel_counts$transcript_id), ]
  novel_frame <- novel_counts[grep("Bambu", novel_counts$transcript_id), ]

  coding_transcripts <- unique(transcriptome_coding[transcriptome_coding$type == "CDS"]$transcript_id)
  sqanti_categories <- vroom::vroom("results/annotate/run_sqanti/kinnex_wtc_11/kinnex_wtc_11_classification.txt") %>%
    filter(grepl("Bambu", isoform)) %>%
    dplyr::select(isoform, structural_category, subcategory) %>%
    rename(isoform = "transcript_id")

  sqanti_categories$structural_category <- ifelse(
    sqanti_categories$structural_category == "novel_in_catalog",
    "NIC",
    ifelse(
      sqanti_categories$structural_category == "novel_not_in_catalog",
      "NNIC",
      ifelse(
        sqanti_categories$structural_category == "fusion",
        "Fusion",
        ifelse(
          sqanti_categories$structural_category == "antisense",
          "Antisense",
          ifelse(
            sqanti_categories$structural_category == "genic",
            "Genic",
            ifelse(
              sqanti_categories$structural_category == "incomplete-splice_match",
              "ISM",
              "Intergenic"
            )
          )
        )
      )
    )
  )


  plt_frame <- data.frame(
    transcript_id = novel_frame$transcript_id,
    cpm = apply(cps, 1, function(x) mean(sort(x, decreasing = TRUE)[1:2])),
    coding = novel_frame$transcript_id %in% coding_transcripts
  ) %>% left_join(sqanti_categories)

  plt_frame$cut_cpm <- cut(plt_frame$cpm, c(-1, 10, max(plt_frame$cpm) + 0.01))

  plt_frame$cut_cpm <- ifelse(
    is.na(plt_frame$cut_cpm),
    levels(plt_frame$cut_cpm)[1],
    levels(plt_frame$cut_cpm)[as.numeric(plt_frame$cut_cpm)]
  )
  plt_frame$cut_cpm <- paste0("CPM: ", plt_frame$cut_cpm)

  plt_frame$cut_cpm <- ifelse(plt_frame$cut_cpm == "CPM: (10,1.56e+03]",
    "CPM > 10",
    "CPM: (0,10]"
  )

  plt_frame$cut_cpm <- factor(
    plt_frame$cut_cpm,
    levels = c(
      "CPM: (0,10]",
      "CPM > 10"
    )
  )

  plt_frame$coding <- ifelse(
    plt_frame$coding,
    "Coding",
    "Non-coding"
  )

  novel_isoform_overview <- plt_frame %>% ggplot(aes(x = structural_category, fill = coding)) +
    geom_bar(position = "stack") +
    facet_wrap(~cut_cpm, nrow = 1, scales = "free_y") +
    ggpubfigs::theme_big_simple() +
    labs(x = "", y = "Novel isoforms", fill = "") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = config$coding_noncoding_cols) +
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
    ) +
    theme(axis.title.x = element_blank())

  fourth_row <- cowplot::plot_grid(
    cowplot::plot_grid(
      cowplot::plot_grid(NULL, cowplot::plot_grid(NULL, mds_gencode, ncol = 2, rel_widths = c(0.05, 0.975)), nrow = 2, rel_heights = c(0.025, 0.975)),
      nrow = 1,
      rel_heights = c(1),
      labels = c("I"),
      label_size = config$label_size
    ),
    cowplot::plot_grid(
      cowplot::plot_grid(NULL,
        cowplot::plot_grid(NULL,
          novel_isoform_overview,
          ncol = 2,
          rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      nrow = 1,
      labels = c("N"),
      label_size = config$label_size
    ),
    cowplot::plot_grid(
      cowplot::plot_grid(NULL,
        cowplot::plot_grid(NULL,
          (ggdraw() + draw_image(magick::image_read_pdf(config$novel_track_path, density = 800))),
          ncol = 2,
          rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      nrow = 1,
      labels = c("O"),
      label_size = config$label_size
    ),
    ncol = 3,
    rel_widths = c(0.45, 0.3, 0.35)
  )

  transcriptome <- import("results/discover/fix_bambu_gene_ids/transcriptome.gtf")

  length_frame <- data.frame(transcriptome) %>%
    filter(type == "exon") %>%
    mutate(exonic_length = end - start) %>%
    group_by(transcript_id) %>%
    summarise(exonic_length = sum(exonic_length))

  kinnex_novel_frame <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method_novel, "_novel/gencode/transcript_counts_formatted.tsv")
  )[, -7]


  s <- edgeR::catchSalmon(
    c(
      "results/quantify/run_salmon_illumina_novel/gencode/day0-rep1",
      "results/quantify/run_salmon_illumina_novel/gencode/day0-rep2",
      "results/quantify/run_salmon_illumina_novel/gencode/day0-rep3",
      "results/quantify/run_salmon_illumina_novel/gencode/day1-rep1",
      "results/quantify/run_salmon_illumina_novel/gencode/day3-rep1",
      "results/quantify/run_salmon_illumina_novel/gencode/day3-rep2",
      "results/quantify/run_salmon_illumina_novel/gencode/day4-rep1",
      "results/quantify/run_salmon_illumina_novel/gencode/day5-rep1",
      "results/quantify/run_salmon_illumina_novel/gencode/day5-rep2",
      "results/quantify/run_salmon_illumina_novel/gencode/day5-rep3"
    )
  )

  illumina_novel_frame <- data.frame(
    transcript_id = rownames(s$counts),
    `day0-rep1` = s$counts[, 1],
    `day0-rep2` = s$counts[, 2],
    `day0-rep3` = s$counts[, 3],
    `day1-rep1` = s$counts[, 4],
    `day3-rep1` = s$counts[, 5],
    `day3-rep2` = s$counts[, 6],
    `day4-rep1` = s$counts[, 7],
    `day5-rep1` = s$counts[, 8],
    `day5-rep2` = s$counts[, 9],
    `day5-rep3` = s$counts[, 10],
    check.names = FALSE
  ) %>%
    left_join(
      data.frame(
        gene_id = kinnex_novel_frame$gene_id,
        transcript_id = kinnex_novel_frame$transcript_id
      ) %>% distinct()
    ) %>%
    dplyr::select(
      gene_id,
      transcript_id,
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
    arrange(desc(transcript_id))

  novel_mask <- grepl("Bambu", kinnex_novel_frame$transcript_id)
  novel_mask_ill <- grepl("Bambu", illumina_novel_frame$transcript_id)
  length_frame_ordered <- length_frame[match(kinnex_novel_frame$transcript_id, length_frame$transcript_id), ]
  n_novel_kinnex <- apply(edgeR::cpm(kinnex_novel_frame[, -(1:2)]), 2, function(x) rep(1, sum(x[which(novel_mask)] > config$cpm_threshold)))
  novel_kinnex_length <- apply(edgeR::cpm(kinnex_novel_frame[, -(1:2)]), 2, function(x) length_frame_ordered$exonic_length[which(x[which(novel_mask)] > config$cpm_threshold)])
  novel_kinnex_cpm <- apply(edgeR::cpm(kinnex_novel_frame[, -(1:2)]), 2, function(x) x[which(novel_mask)][which(x[which(novel_mask)] > config$cpm_threshold)])

  n_novel_illumina <- apply(edgeR::cpm(illumina_novel_frame[, -(1:2)]), 2, function(x) rep(1, sum(x[which(novel_mask_ill)] > config$cpm_threshold)))
  novel_illumina_length <- apply(edgeR::cpm(illumina_novel_frame[, -(1:2)]), 2, function(x) length_frame_ordered$exonic_length[which(x[which(novel_mask_ill)] > config$cpm_threshold)])
  novel_illumina_cpm <- apply(edgeR::cpm(illumina_novel_frame[, -(1:2)]), 2, function(x) x[which(novel_mask)][which(x[which(novel_mask_ill)] > config$cpm_threshold)])

  n_kinnex <- apply(edgeR::cpm(kinnex_novel_frame[, -(1:2)]), 2, function(x) rep(1, sum(x[which(!novel_mask)] > config$cpm_threshold)))
  kinnex_length <- apply(edgeR::cpm(kinnex_novel_frame[, -(1:2)]), 2, function(x) length_frame_ordered$exonic_length[which(x[which(!novel_mask)] > config$cpm_threshold)])
  kinnex_cpm <- apply(edgeR::cpm(kinnex_novel_frame[, -(1:2)]), 2, function(x) x[which(!novel_mask)][which(x[which(!novel_mask)] > config$cpm_threshold)])

  n_illumina <- apply(edgeR::cpm(illumina_novel_frame[, -(1:2)]), 2, function(x) rep(1, sum(x[which(!novel_mask_ill)] > config$cpm_threshold)))
  illumina_length <- apply(edgeR::cpm(illumina_novel_frame[, -(1:2)]), 2, function(x) length_frame_ordered$exonic_length[which(x[which(!novel_mask_ill)] > config$cpm_threshold)])
  illumina_cpm <- apply(edgeR::cpm(illumina_novel_frame[, -(1:2)]), 2, function(x) x[which(!novel_mask_ill)][which(x[which(!novel_mask_ill)] > config$cpm_threshold)])


  plt_frame <- data.frame(
    n = c(
      unlist(n_novel_kinnex),
      unlist(n_novel_illumina),
      unlist(n_kinnex),
      unlist(n_illumina)
    ),
    length = c(
      unlist(novel_kinnex_length),
      unlist(novel_illumina_length),
      unlist(kinnex_length),
      unlist(illumina_length)
    ),
    tech = c(
      rep("Kinnex", length(unlist(n_novel_kinnex))),
      rep("Illumina", length(unlist(n_novel_illumina))),
      rep("Kinnex", length(unlist(n_kinnex))),
      rep("Illumina", length(unlist(n_illumina)))
    ),
    novel = c(
      rep("Novel (from Kinnex)", length(unlist(n_novel_kinnex))),
      rep("Novel (from Kinnex)", length(unlist(n_novel_illumina))),
      rep("Reference", length(unlist(n_kinnex))),
      rep("Reference", length(unlist(n_illumina)))
    ),
    day = c(
      rep(names(n_novel_kinnex), sapply(n_novel_kinnex, length)),
      rep(names(n_novel_illumina), sapply(n_novel_illumina, length)),
      rep(names(n_kinnex), sapply(n_kinnex, length)),
      rep(names(n_illumina), sapply(n_illumina, length))
    ),
    cpm = c(
      unlist(novel_kinnex_cpm),
      unlist(novel_illumina_cpm),
      unlist(kinnex_cpm),
      unlist(illumina_cpm)
    )
  )

  plt_frame$cut_length <- (cut(
    plt_frame$length,
    c(7, 500, 1000, 1000 * 10, 347559)
  ))

  plt_frame$cut_cpm <- ((cut(
    plt_frame$cpm,
    c(-1, 50, 21531.82 + 0.01)
  )))


  plt_frame$facet <- paste0(
    plt_frame$tech, ifelse(
      plt_frame$cut_cpm == "(-1,50]",
      " (CPM <= 50)",
      " (CPM > 50)"
    )
  )

  expressed_transcripts_low <- plt_frame %>%
    filter(cut_cpm == "(-1,50]") %>%
    ggplot(aes(fill = novel, x = day)) +
    geom_bar(position = "stack") +
    facet_wrap(~facet, nrow = 1) +
    theme_big_simple() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = "", x = "Sample", y = "Transcripts") +
    scale_fill_manual(values = rev(config$known_novel_cols))

  expressed_transcripts_high <- plt_frame %>%
    filter(cut_cpm != "(-1,50]") %>%
    ggplot(aes(fill = novel, x = day)) +
    geom_bar(position = "stack") +
    facet_wrap(~facet, nrow = 1) +
    theme_big_simple() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = "", x = "Sample", y = "Transcripts") +
    scale_fill_manual(values = rev(config$known_novel_cols))

  kinnex_novel_frame <- vroom::vroom(
    paste0("results/format_quantify/", config$chosen_method_novel, "_novel/gencode/gene_counts_formatted.tsv")
  )[, -6]

  illumina_novel_frame <- vroom::vroom(
    "results/format_quantify/salmon_illumina_novel/gencode/gene_counts_formatted.tsv"
  )[, -6]


  txdb <- makeTxDbFromGFF("results/discover/fix_bambu_gene_ids/transcriptome.gtf", format = "gtf")
  exons_per_gene <- exonsBy(txdb, by = "gene")
  summed_exonic_gene_size <- sum(width(reduce(exons_per_gene)))
  length_frame <- data.frame(
    gene_id = names(summed_exonic_gene_size),
    exonic_length = unname(summed_exonic_gene_size)
  )

  novel_mask <- grepl("Bambu", kinnex_novel_frame$gene_id)

  length_frame_ordered <- length_frame[match(kinnex_novel_frame$gene_id, length_frame$gene_id), ]

  n_novel_kinnex <- apply(edgeR::cpm(kinnex_novel_frame[, -(1)]), 2, function(x) rep(1, sum(x[which(novel_mask)] > config$cpm_threshold)))
  novel_kinnex_length <- apply(edgeR::cpm(kinnex_novel_frame[, -(1)]), 2, function(x) length_frame_ordered$exonic_length[which(x[which(novel_mask)] > config$cpm_threshold)])
  novel_kinnex_cpm <- apply(edgeR::cpm(kinnex_novel_frame[, -(1)]), 2, function(x) x[which(novel_mask)][which(x[which(novel_mask)] > config$cpm_threshold)])

  n_novel_illumina <- apply(edgeR::cpm(illumina_novel_frame[, -(1)]), 2, function(x) rep(1, sum(x[which(novel_mask)] > config$cpm_threshold)))
  novel_illumina_length <- apply(edgeR::cpm(illumina_novel_frame[, -(1)]), 2, function(x) length_frame_ordered$exonic_length[which(x[which(novel_mask)] > config$cpm_threshold)])
  novel_illumina_cpm <- apply(edgeR::cpm(illumina_novel_frame[, -(1)]), 2, function(x) x[which(novel_mask)][which(x[which(novel_mask)] > config$cpm_threshold)])

  n_kinnex <- apply(edgeR::cpm(kinnex_novel_frame[, -(1)]), 2, function(x) rep(1, sum(x[which(!novel_mask)] > config$cpm_threshold)))
  kinnex_length <- apply(edgeR::cpm(kinnex_novel_frame[, -(1)]), 2, function(x) length_frame_ordered$exonic_length[which(x[which(!novel_mask)] > config$cpm_threshold)])
  kinnex_cpm <- apply(edgeR::cpm(kinnex_novel_frame[, -(1)]), 2, function(x) x[which(!novel_mask)][which(x[which(!novel_mask)] > config$cpm_threshold)])

  n_illumina <- apply(edgeR::cpm(illumina_novel_frame[, -(1)]), 2, function(x) rep(1, sum(x[which(!novel_mask)] > config$cpm_threshold)))
  illumina_length <- apply(edgeR::cpm(illumina_novel_frame[, -(1)]), 2, function(x) length_frame_ordered$exonic_length[which(x[which(!novel_mask)] > config$cpm_threshold)])
  illumina_cpm <- apply(edgeR::cpm(illumina_novel_frame[, -(1)]), 2, function(x) x[which(!novel_mask)][which(x[which(!novel_mask)] > config$cpm_threshold)])

  plt_frame <- data.frame(
    n = c(
      unlist(n_novel_kinnex),
      unlist(n_novel_illumina),
      unlist(n_kinnex),
      unlist(n_illumina)
    ),
    length = c(
      unlist(novel_kinnex_length),
      unlist(novel_illumina_length),
      unlist(kinnex_length),
      unlist(illumina_length)
    ),
    tech = c(
      rep("Kinnex", length(unlist(n_novel_kinnex))),
      rep("Illumina", length(unlist(n_novel_illumina))),
      rep("Kinnex", length(unlist(n_kinnex))),
      rep("Illumina", length(unlist(n_illumina)))
    ),
    novel = c(
      rep("Novel (from Kinnex)", length(unlist(n_novel_kinnex))),
      rep("Novel (from Kinnex)", length(unlist(n_novel_illumina))),
      rep("Reference", length(unlist(n_kinnex))),
      rep("Reference", length(unlist(n_illumina)))
    ),
    day = c(
      rep(names(n_novel_kinnex), sapply(n_novel_kinnex, length)),
      rep(names(n_novel_illumina), sapply(n_novel_illumina, length)),
      rep(names(n_kinnex), sapply(n_kinnex, length)),
      rep(names(n_illumina), sapply(n_illumina, length))
    ),
    cpm = c(
      unlist(novel_kinnex_cpm),
      unlist(novel_illumina_cpm),
      unlist(kinnex_cpm),
      unlist(illumina_cpm)
    )
  )

  plt_frame %>%
    filter(tech == "Illumina") %>%
    group_by(novel) %>%
    summarise(mean = n() / 10)

  plt_frame$cut_length <- (cut(
    plt_frame$length,
    c(7, 500, 1000, 1000 * 10, 347559)
  ))

  quantile(plt_frame$cpm)

  plt_frame$cut_cpm <- ((cut(
    plt_frame$cpm,
    c(-1, 50, 21610.724939 + 0.01)
  )))

  plt_frame$facet <- paste0(
    plt_frame$tech, ifelse(
      plt_frame$cut_cpm == "(-1,50]",
      " (CPM <= 50)",
      " (CPM > 50)"
    )
  )

  expressed_genes_low <- plt_frame %>%
    filter(cut_cpm == "(-1,50]") %>%
    ggplot(aes(fill = novel, x = day)) +
    geom_bar(position = "stack") +
    facet_wrap(~facet, nrow = 1) +
    theme_big_simple() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = "", x = "Sample", y = "Genes") +
    scale_fill_manual(values = rev(config$known_novel_cols))

  expressed_genes_high <- plt_frame %>%
    filter(cut_cpm != "(-1,50]") %>%
    ggplot(aes(fill = novel, x = day)) +
    geom_bar(position = "stack") +
    facet_wrap(~facet, nrow = 1) +
    theme_big_simple() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = "", x = "Sample", y = "Genes") +
    scale_fill_manual(values = rev(config$known_novel_cols))

  fifth_row <- cowplot::plot_grid(
    cowplot::plot_grid(
      cowplot::plot_grid(NULL, cowplot::plot_grid(NULL,
        vol_plot,
        ncol = 2,
        rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      cowplot::plot_grid(NULL, cowplot::plot_grid(NULL,
        spline_plot,
        ncol = 2,
        rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
      ),
      nrow = 2,
      label_size = config$label_size,
      labels = c("L", "M")
    ),
    cowplot::plot_grid(
      cowplot::plot_grid(
        cowplot::plot_grid(NULL, cowplot::plot_grid(NULL,
          expressed_genes_low + theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          ) + theme(legend.position = "none"),
          ncol = 2,
          rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
        ),
        cowplot::plot_grid(NULL, cowplot::plot_grid(NULL,
          expressed_genes_high + theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          ) + theme(legend.position = "none"),
          ncol = 2,
          rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
        ),
        ncol = 2
      ),
      cowplot::plot_grid(
        cowplot::plot_grid(NULL, cowplot::plot_grid(NULL,
          expressed_transcripts_low,
          ncol = 2,
          rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
        ),
        cowplot::plot_grid(NULL, cowplot::plot_grid(NULL,
          expressed_transcripts_high,
          ncol = 2,
          rel_widths = c(config$side_difference_large, (1 - config$side_difference_large))
        ),
        nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
        ),
        ncol = 2
      ),
      nrow = 2,
      rel_heights = c(0.4, 0.75)
    ),
    label_size = config$label_size,
    ncol = 2,
    labels = c("", "P"),
    rel_widths = c(0.425, 0.6)
  )

  final <- cowplot::plot_grid(
    first_row,
    second_row,
    third_row,
    fourth_row,
    fifth_row,
    nrow = 5,
    rel_heights = c(1, 0.8, 1, 1 / 2, 1.1)
  )

  ggsave(output_path, final,
    dpi = 300,
    width = 26,
    height = 36
  )
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(data.table)
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
  library(RColorBrewer)
  library(ggplotify)
  library(rjson)
  library(GenomicFeatures)
})

sessionInfo()
status <- plot_figure_01(
  output_path = snakemake@output[[1]]
)

sink()
sink()

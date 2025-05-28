plot_figure_S03 <- function(output_path, depth) {
  config <- rjson::fromJSON(file = "config/config.json")
  params <- rjson::fromJSON(file = "config/params.json")

  plt_frame <- data.frame(
    tech = c(
      rep("Kinnex", 8),
      rep("Illumina", 8)
    ),
    metric_type = c(
      rep(rep(c("Recall", "Precision"), each = 4), 2)
    ),
    value = c(
      params$discover_sirv_recall_kinnex,
      params$discover_sirv_precision_kinnex,
      params$discover_sirv_recall_illumina,
      params$discover_sirv_precision_illumina
    ),
    sizes = rep(paste0(as.character(params$discover_sirv_sizes / 1e6), " M"), 4)
  )

  plt_frame$sizes <- paste0("Reads: ", plt_frame$sizes)

  plt_frame$sizes <- factor(
    plt_frame$sizes,
    levels = c(
      "Reads: 0.25 M",
      "Reads: 0.5 M",
      "Reads: 1 M",
      "Reads: 2.5 M"
    )
  )

  discovery_performance <- plt_frame %>%
    filter(sizes == "Reads: 2.5 M") %>%
    ggplot(aes(x = metric_type, y = value, fill = tech)) +
    geom_bar(position = "dodge", stat = "identity") +
    ggpubfigs::theme_big_simple() +
    labs(x = "", y = "Performance", fill = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = rev(c("#D8178C", "#ffb441"))) +
    theme(axis.title.x = element_blank())

  kinnex_detected <- vroom::vroom("results/discover_sirv/evaluate_performance/pb_1_2500000.0.tracking", col_names = FALSE, delim = "\t") %>% filter(X4 == "=")
  illumina_detected <- vroom::vroom("results/discover_sirv/evaluate_performance/illumina_1_2500000.0.tracking", col_names = FALSE, delim = "\t") %>% filter(X4 == "=")

  union_detected <- union(
    sapply(strsplit(kinnex_detected$X3, "\\|"), function(x) x[[2]]),
    sapply(strsplit(illumina_detected$X3, "\\|"), function(x) x[[2]])
  )

  coldata_ill <- data.frame(
    files = c(
      "results/quantify/run_salmon_illumina/sirv/day0-rep1/quant.sf",
      "results/quantify/run_salmon_illumina/sirv/day0-rep2/quant.sf",
      "results/quantify/run_salmon_illumina/sirv/day0-rep3/quant.sf",
      "results/quantify/run_salmon_illumina/sirv/day5-rep1/quant.sf",
      "results/quantify/run_salmon_illumina/sirv/day5-rep2/quant.sf",
      "results/quantify/run_salmon_illumina/sirv/day5-rep3/quant.sf"
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
  illumina_se <- fishpond::computeInfRV(se)[93:161, ]

  plt_frame <- data.frame(
    transcript_id = rownames(illumina_se),
    inf_var = rowData(illumina_se)$meanInfRV,
    discovered_in = ifelse(
      rownames(illumina_se) %in% sapply(strsplit(kinnex_detected$X3, "\\|"), function(x) x[[2]]) & rownames(illumina_se) %in% sapply(strsplit(illumina_detected$X3, "\\|"), function(x) x[[2]]),
      "Both",
      ifelse(
        rownames(illumina_se) %in% sapply(strsplit(kinnex_detected$X3, "\\|"), function(x) x[[2]]),
        "Kinnex",
        ifelse(
          rownames(illumina_se) %in% sapply(strsplit(illumina_detected$X3, "\\|"), function(x) x[[2]]),
          "Illumina",
          "Neither"
        )
      )
    )
  )

  set.seed(config$seed)

  discovery_performance_inf_var <- plt_frame %>% ggplot(aes(x = discovered_in, y = inf_var, color = discovered_in)) +
    geom_jitter(size = 2.5) +
    theme_big_simple() +
    labs(
      color = "Discovered by", y = "Inferential variability"
    ) +
    scale_color_manual(values = c("red", "#ffb441", "#D9188D", "grey")) +
    theme(axis.title.x = element_blank())

  all_counts <- Reduce(c, lapply(
    c("oarfish"), function(method) {
      lapply(
        c("5000000.0"), function(depth) {
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
    c(config$chosen_method, "salmon_illumina"), function(method) {
      lapply(
        c("30000000.0"), function(depth) {
          vroom::vroom(
            paste0("results/format_quantify_subsampled/", method, "/gencode/1/", depth, "/transcript_counts_formatted.tsv")
          )
        }
      )
    }
  ))

  dte_calls <- lapply(
    all_counts,
    function(x) {
      dge <- edgeR::DGEList(counts = x[, -(1:2)], genes = x[, 1:2], group = c(0, 0, 0, 1, 1, 1))


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

      return(data.frame(
        transcript_id = data.frame(tt) %>% arrange(desc(transcript_id)) %>% pull(transcript_id),
        fdr = data.frame(tt) %>% arrange(desc(transcript_id)) %>% pull(FDR)
      ))
    }
  )

  joint_transcripts <- intersect(dte_calls[[1]]$transcript_id, dte_calls[[2]]$transcript_id)
  plt_frame <- data.frame(
    transcript_id = joint_transcripts,
    kinnex = dte_calls[[1]]$fdr[match(joint_transcripts, dte_calls[[1]]$transcript_id)],
    ill = dte_calls[[2]]$fdr[match(joint_transcripts, dte_calls[[2]]$transcript_id)],
    correction = "No correction"
  )

  all_counts <- Reduce(c, lapply(
    c(config$chosen_method, "run_salmon_illumina_bootstrap"), function(method) {
      lapply(
        c("30000000.0"), function(depth) {
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

  dte_calls <- lapply(
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

      return(data.frame(
        transcript_id = data.frame(tt) %>% arrange(desc(transcript_id)) %>% pull(transcript_id),
        fdr = data.frame(tt) %>% arrange(desc(transcript_id)) %>% pull(FDR)
      ))
    }
  )

  joint_transcripts <- intersect(dte_calls[[1]]$transcript_id, dte_calls[[2]]$transcript_id)
  plt_frame <- rbind(
    plt_frame,
    data.frame(
      transcript_id = joint_transcripts,
      kinnex = dte_calls[[1]]$fdr[match(joint_transcripts, dte_calls[[1]]$transcript_id)],
      ill = dte_calls[[2]]$fdr[match(joint_transcripts, dte_calls[[2]]$transcript_id)],
      correction = "With correction"
    )
  )

  plt_frame <- plt_frame %>% filter(transcript_id %in% names(which(table(plt_frame$transcript_id) == 2)))

  plt_frame$correction <- factor(
    plt_frame$correction,
    levels = c("No correction", "With correction")
  )



  plt_frame$sig <- ifelse(
    plt_frame$ill < 0.01 & plt_frame$kinnex < 0.01,
    "Both",
    ifelse(
      plt_frame$ill < 0.01,
      "Illumina",
      ifelse(
        plt_frame$kinnex < 0.01,
        "Kinnex",
        "Neither"
      )
    )
  )

  fdr_conc_dte <- plt_frame %>% ggplot(aes(x = -log(ill, base = 10), y = -log(kinnex, base = 10), color = sig)) +
    geom_point(size = 2.5) +
    facet_wrap(~correction) +
    geom_hline(yintercept = 2, lty = 2, linewidth = 1) +
    geom_vline(xintercept = 2, lty = 2, linewidth = 1) +
    theme_big_simple() +
    scale_color_manual(values = c("red", "#ffb441", "#D9188D", "grey")) +
    labs(x = "-log10(FDR) (Illumina)", y = "-log10(FDR) (Kinnex)", color = "Significance") +
    ggtitle("GENCODE DTE: FDR concordance with and without correction")

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
    filter(Illumina < 1) %>%
    pull(transcript_id)

  all_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("salmon_illumina", "bambu", "isoquant", "oarfish", "kallisto", "salmon"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/gencode/", rep, "/30000000.0/transcript_counts_formatted.tsv")
        ) %>%
          filter(transcript_id %in% non_variable_transcripts) %>%
          arrange(desc(transcript_id))
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
    which(apply(prepped_counts[[1]], 1, function(x) sum(x[1:3] > 1.0) >= 3)),
    which(apply(prepped_counts[[4]], 1, function(x) sum(x[1:3] > 1.0) >= 3))
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
      5
    ),
    rep = rep(rep(1:3, each = length(mask)), 5),
    method = rep(
      c(
        "Kinnex (Bambu)",
        "Kinnex (Isoquant)",
        "Kinnex (Oarfish)",
        "Kinnex (Kallisto)",
        "Kinnex (Salmon)"
      ),
      each = length(mask) * 3
    )
  )

  plt_frame$length_factor <- factor(cut(plt_frame$length, c(0, 1250, 5000, 347559 + 1)))

  plt_frame$illumina_expr_factor <- factor(cut(plt_frame$illumina_fc, c(-9, -0.7412485, -0.1025366, 0.4798768, 8.3373789 + 1)))


  plt_frame$length_factor <- paste0(
    "Length: ", plt_frame$length_factor
  )

  plt_frame$length_factor <- ifelse(
    plt_frame$length_factor == "Length: (5e+03,3.48e+05]",
    "Length > 3.48e+05",
    plt_frame$length_factor
  )

  plt_frame$length_factor <- factor(
    plt_frame$length_factor,
    levels = c(
      "Length: (0,1.25e+03]",
      "Length: (1.25e+03,5e+03]",
      "Length > 3.48e+05"
    )
  )
  set.seed(42)

  illumina_relative <- plt_frame %>%
    filter(method == "Kinnex (Oarfish)") %>%
    ggplot(aes(x = illumina_fc, y = log_fc, shape = as.factor(rep))) +
    geom_point(size = 4, alpha = .7) +
    labs(x = "logFC (Illumina)", y = "logFC (Kinnex)", shape = "Replicate") +
    ggpubfigs::theme_big_simple() +
    stat_cor(aes(label = ..r.label.., colour = NULL, group = 1), show.legend = FALSE, method = "spearman", size = 6, label.x.npc = c(0.025), label.y.npc = c(1.0)) +
    facet_wrap(~length_factor) +
    geom_density_2d(color = "lightblue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

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


  transcript_frame <- transcript_frame %>% left_join(
    data.frame(
      transcript_id = names(rowData(kinnex_se)$meanInfRV),
      inf_var = unname(rowData(kinnex_se)$meanInfRV)
    )
  )

  transcriptome_strings <- readDNAStringSet("results/prepare/extract_transcriptomes/gencode_transcriptome.fa")

  gc_content <- letterFrequency(transcriptome_strings, letters = "GC", as.prob = TRUE)

  transcript_frame <- transcript_frame %>% left_join(
    data.frame(
      transcript_id = names(transcriptome_strings),
      gc = gc_content[, 1]
    )
  )


  plt_frame <- data.frame(
    transcript_id = transcript_frame$transcript_id,
    inf_var = transcript_frame$inf_var,
    gc = transcript_frame$gc,
    length = transcript_frame$length,
    observed_abundance = unlist(lapply(
      prepped_counts[-1],
      function(counts) {
        c(
          counts[, 1],
          counts[, 2],
          counts[, 3],
          counts[, 4],
          counts[, 5],
          counts[, 6]
        )
      }
    )),
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

  plt_frame$length_factor <- factor(cut(plt_frame$length, c(0, 1250, 5000, 347559 + 1)))


  plt_frame$length_factor <- paste0(
    "Length: ", plt_frame$length_factor
  )

  plt_frame$length_factor <- ifelse(
    plt_frame$length_factor == "Length: (5e+03,3.48e+05]",
    "Length > 3.48e+05",
    plt_frame$length_factor
  )

  plt_frame$length_factor <- factor(
    plt_frame$length_factor,
    levels = c(
      "Length: (0,1.25e+03]",
      "Length: (1.25e+03,5e+03]",
      "Length > 3.48e+05"
    )
  )

  plt_frame$inf_var_factor <- factor(cut(plt_frame$inf_var, c(0, 0.1, 165)))

  plt_frame$inf_var_factor <- paste0(
    "Inf. variability: ", plt_frame$inf_var_factor
  )

  illumina_absolute <- plt_frame %>%
    filter(condition == "Day 0" & method == "Kinnex (Oarfish)") %>%
    ggplot(aes(x = illumina_abundance, y = observed_abundance, shape = as.factor(rep))) +
    geom_point(size = 4, alpha = .7) +
    geom_density_2d(color = "lightblue") +
    labs(x = "logCPM (Illumina)", y = "logCPM (Kinnex)", shape = "Replicate") +
    ggpubfigs::theme_big_simple() +
    stat_cor(aes(label = ..r.label.., colour = NULL, group = 1), show.legend = FALSE, method = "spearman", size = 6, label.x.npc = c(0.025), label.y.npc = c(1)) +
    facet_wrap(~length_factor) +
    scale_color_manual(values = ggpubfigs::friendly_pals$zesty_four[c(1, 3)]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  
  
  
  oarfish_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("oarfish"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/gencode/", rep, "/30000000.0/transcript_counts_formatted.tsv")
        ) %>%
          arrange(desc(transcript_id))
      })
    }
  ))[[1]]
  
  illumina_counts <- Reduce(c, lapply(
    1:1,
    function(rep) {
      lapply(c("salmon_illumina"), function(method) {
        vroom::vroom(
          paste0("results/format_quantify_subsampled/", method, "/gencode/", rep, "/30000000.0/transcript_counts_formatted.tsv")
        ) %>%
          arrange(desc(transcript_id))
      })
    }
  ))[[1]]
  
  coldata_ill <- data.frame(
    files = c(
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/30000000.0/day0-rep1/quant.sf"),
      paste0("results/quantify_subsampled/run_salmon_illumina_bootstrap/gencode/1/30000000.0/day0-rep2/quant.sf"),
      paste0("results/quantify_subsampled//run_salmon_illumina_bootstrap/gencode/1/30000000.0/day0-rep3/quant.sf")
    ),
    names = c(
      "day0-rep1",
      "day0-rep2",
      "day0-rep3"
    )
  )
  
  se <- tximeta(coldata_ill)
  illumina_se <- fishpond::computeInfRV(se)
  
  
  
  
  
  
  plt_frame <- data.frame(
    tech = c(rep("Kinnex", 6), rep("Illumina", 6 * 3)),
    transcript_id = rep(
      rep(
        c(
          "ENST00000455834",
          "ENST00000448997"
        ),
        each = 3
      ), 4
    ),
    cpm = c(
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000114735.10", ]$transcript_id), ][2, 1:3],
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000177889.10", ]$transcript_id), ][3, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000455834.5",
        "ENST00000448997.5"
      )), ][2, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000455834.5",
        "ENST00000448997.5"
      )), ][1, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000455834.5",
        "ENST00000448997.5"
      )), ][2, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000455834.5",
        "ENST00000448997.5"
      )), ][1, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000455834.5",
        "ENST00000448997.5"
      )), ][2, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000455834.5",
        "ENST00000448997.5"
      )), ][1, 1:3]
    ),
    rep = c(
      rep(c("Day 0-1", "Day 0-2", "Day 0-3"), 4),
      rep(c("Day 0-1 (Inf. Rep. 1)", "Day 0-2 (Inf. Rep. 1)", "Day 0-3 (Inf. Rep. 1)"), 2),
      rep(c("Day 0-1 (Inf. Rep. 2)", "Day 0-2 (Inf. Rep. 2)", "Day 0-3 (Inf. Rep. 2)"), 2)
    )
  )
  
  table(plt_frame$rep)
  
  plt_frame$rep <- factor(
    plt_frame$rep,
    levels = c(
      "Day 0-1",
      "Day 0-2",
      "Day 0-3",
      "Day 0-1 (Inf. Rep. 1)",
      "Day 0-2 (Inf. Rep. 1)",
      "Day 0-3 (Inf. Rep. 1)",
      "Day 0-1 (Inf. Rep. 2)",
      "Day 0-2 (Inf. Rep. 2)",
      "Day 0-3 (Inf. Rep. 2)"
    )
  )
  
  supp_example_01 <- plt_frame %>% ggplot(aes(x = transcript_id, fill = tech, y = cpm + 1)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~rep, nrow = 1) +
    theme_big_simple() +
    scale_fill_manual(
      values = c("#ffb441", "#D8178C")
    ) +
    scale_y_log10() +
    labs(y = "CPM", x = "", fill = "", title = "HEMK1 - Inf. var.: 0.39 -log10(q): 5.7 CPM: 6.47") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_blank())
  
  plt_frame <- data.frame(
    tech = c(rep("Kinnex", 6), rep("Illumina", 6 * 3)),
    transcript_id = rep(
      rep(
        c(
          "ENST00000372004",
          "ENST00000371998"
        ),
        each = 3
      ), 4
    ),
    cpm = c(
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000124151.19", ]$transcript_id), ][3, 1:3],
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000124151.19", ]$transcript_id), ][4, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000372004.7",
        "ENST00000371998.8"
      )), ][2, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000372004.7",
        "ENST00000371998.8"
      )), ][1, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000372004.7",
        "ENST00000371998.8"
      )), ][2, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000372004.7",
        "ENST00000371998.8"
      )), ][1, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000372004.7",
        "ENST00000371998.8"
      )), ][2, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000372004.7",
        "ENST00000371998.8"
      )), ][1, 1:3]
    ),
    rep = c(
      rep(c("Day 0-1", "Day 0-2", "Day 0-3"), 4),
      rep(c("Day 0-1 (Inf. Rep. 1)", "Day 0-2 (Inf. Rep. 1)", "Day 0-3 (Inf. Rep. 1)"), 2),
      rep(c("Day 0-1 (Inf. Rep. 2)", "Day 0-2 (Inf. Rep. 2)", "Day 0-3 (Inf. Rep. 2)"), 2)
    )
  )
  
  table(plt_frame$rep)
  
  plt_frame$rep <- factor(
    plt_frame$rep,
    levels = c(
      "Day 0-1",
      "Day 0-2",
      "Day 0-3",
      "Day 0-1 (Inf. Rep. 1)",
      "Day 0-2 (Inf. Rep. 1)",
      "Day 0-3 (Inf. Rep. 1)",
      "Day 0-1 (Inf. Rep. 2)",
      "Day 0-2 (Inf. Rep. 2)",
      "Day 0-3 (Inf. Rep. 2)"
    )
  )
  
  supp_example_02 <- plt_frame %>% ggplot(aes(x = transcript_id, fill = tech, y = cpm + 1)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~rep, nrow = 1) +
    theme_big_simple() +
    scale_fill_manual(
      values = c("#ffb441", "#D8178C")
    ) +
    scale_y_log10() +
    labs(y = "CPM", x = "", fill = "", title = "NCOA3 - Inf. var.: 52.0 -log10(q): 0.06 CPM: 47.64") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_blank())
  
  plt_frame <- data.frame(
    tech = c(rep("Kinnex", 5 * 3), rep("Illumina", 5 * 3 * 3)),
    transcript_id = rep(
      rep(
        c(
          "ENST00000304567",
          "ENST00000360566",
          "ENST00000485717",
          "ENST00000641498",
          "ENST00000652660"
        ),
        each = 3
      ), 4
    ),
    cpm = c(
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000171848.16", ]$transcript_id), ][17, 1:3],
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000171848.16", ]$transcript_id), ][16, 1:3],
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000171848.16", ]$transcript_id), ][11, 1:3],
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000171848.16", ]$transcript_id), ][4, 1:3],
      edgeR::cpm(oarfish_counts[, -(1:2)])[which(oarfish_counts$transcript_id %in% oarfish_counts[oarfish_counts$gene_id == "ENSG00000171848.16", ]$transcript_id), ][1, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][3, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][1, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][5, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][4, 1:3],
      edgeR::cpm(assays(illumina_se)[[1]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][2, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][3, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][1, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][5, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][4, 1:3],
      edgeR::cpm(assays(illumina_se)[[4]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][2, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][3, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][1, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][5, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][4, 1:3],
      edgeR::cpm(assays(illumina_se)[[5]])[which(rownames(illumina_se) %in% c(
        "ENST00000304567.10",
        "ENST00000360566.6",
        "ENST00000485717.1",
        "ENST00000641498.1",
        "ENST00000652660.1"
      )), ][2, 1:3]
    ),
    rep = c(
      rep(c("Day 0-1", "Day 0-2", "Day 0-3"), 5 * 2),
      rep(c("Day 0-1 (Inf. Rep. 1)", "Day 0-2 (Inf. Rep. 1)", "Day 0-3 (Inf. Rep. 1)"), 5),
      rep(c("Day 0-1 (Inf. Rep. 2)", "Day 0-2 (Inf. Rep. 2)", "Day 0-3 (Inf. Rep. 2)"), 5)
    )
  )
  table(plt_frame$rep)
  
  plt_frame$rep <- factor(
    plt_frame$rep,
    levels = c(
      "Day 0-1",
      "Day 0-2",
      "Day 0-3",
      "Day 0-1 (Inf. Rep. 1)",
      "Day 0-2 (Inf. Rep. 1)",
      "Day 0-3 (Inf. Rep. 1)",
      "Day 0-1 (Inf. Rep. 2)",
      "Day 0-2 (Inf. Rep. 2)",
      "Day 0-3 (Inf. Rep. 2)"
    )
  )
  
  supp_example_03 <- plt_frame %>% ggplot(aes(x = transcript_id, fill = tech, y = cpm + 1)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~rep, nrow = 1) +
    theme_big_simple() +
    scale_fill_manual(
      values = c("#ffb441", "#D8178C")
    ) +
    scale_y_log10() +
    labs(y = "CPM", x = "", fill = "", title = "RRM2 - Inf. var.: 150.5 -log10(q): 7.35 CPM: 263.6551") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_blank())
  

  first_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        discovery_performance,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large + 0.025, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        discovery_performance_inf_var,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large + 0.025, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        fdr_conc_dte,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 3,
    labels = c("A", "B", "C"),
    label_size = 36,
    rel_widths = c(1, 1, 2)
  )


  second_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        illumina_relative,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    cowplot::plot_grid(NULL,
      cowplot::plot_grid(NULL,
        illumina_absolute,
        ncol = 2,
        rel_widths = c(
          config$side_difference_large, (1 - config$side_difference_large)
        )
      ),
      nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 2,
    labels = c("D", "E"),
    label_size = 36,
    rel_widths = c(1, 1)
  )
  
  third_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
                       cowplot::plot_grid(NULL,
                                          supp_example_01,
                                          ncol = 2,
                                          rel_widths = c(
                                            config$side_difference / 4, (1 - config$side_difference)
                                          )
                       ),
                       nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 1,
    labels = c("F"),
    label_size = 36,
    rel_widths = c(1)
  )
  
  fourth_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
                       cowplot::plot_grid(NULL,
                                          supp_example_02,
                                          ncol = 2,
                                          rel_widths = c(
                                            config$side_difference / 4, (1 - config$side_difference)
                                          )
                       ),
                       nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 1,
    labels = c("G"),
    label_size = 36,
    rel_widths = c(1)
  )
  
  fifth_row <- cowplot::plot_grid(
    cowplot::plot_grid(NULL,
                       cowplot::plot_grid(NULL,
                                          supp_example_03,
                                          ncol = 2,
                                          rel_widths = c(
                                            config$side_difference / 4, (1 - config$side_difference)
                                          )
                       ),
                       nrow = 2, rel_heights = c(config$top_difference, (1 - config$top_difference))
    ),
    ncol = 1,
    labels = c("H"),
    label_size = 36,
    rel_widths = c(1)
  )

  final <- cowplot::plot_grid(
    first_row,
    second_row,
    third_row,
    fourth_row,
    fifth_row,
    nrow = 5,
    rel_heights = c(1, 1, 1, 1, 1)
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
  library(fishpond)
  library(tximeta)
  library(tximport)
  library(ComplexUpset)
  library(magrittr)
  library(tibble)
  library(Biostrings)
})

sessionInfo()
status <- plot_figure_S03(
  output_path = snakemake@output[[1]],
  depth = snakemake@params[["depth"]]
)

sink()
sink()

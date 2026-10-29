if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# 输入和输出目录 --------------------------------------------------------------
argv <- commandArgs(TRUE)
if(length(argv)!=2) stop("Usage: Rscript SCRIPT ANALYSIS_OUTPUT_DIR NEW_FIGURE_DIR")
input_dir <- argv[1]
output_dir <- argv[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rd_file <- file.path(input_dir, "metaDataDF.csv")
cnv_file <- file.path(input_dir, "CNVresult_filter.csv")
scaffold_file <- file.path(input_dir, "CNVscaffold_background.csv")

sample_names <- c("FASR", "MCF7", "TAMR")
autosomes <- paste0("chr", 1:22)

# 按基因组坐标定义窗口。两个版本共用100 Kb窗口，仅步长不同。
window_bp <- 100000
step_values <- c(100000)
show_unassessed <- FALSE

required_files <- c(rd_file, cnv_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("缺少输入文件: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

# 读取并检查RD数据 ------------------------------------------------------------
rd <- fread(rd_file, encoding = "UTF-8", showProgress = TRUE)
required_rd_cols <- c("Chromosome", "BinLowEdge", "BinUpEdge", sample_names)
missing_rd_cols <- setdiff(required_rd_cols, names(rd))
if (length(missing_rd_cols) > 0L) {
  stop("metaDataDF.csv缺少列: ", paste(missing_rd_cols, collapse = ", "),
       call. = FALSE)
}

rd[, source_row := .I]
rd <- rd[Chromosome %chin% autosomes]
rd[, Chromosome := factor(as.character(Chromosome), levels = autosomes, ordered = TRUE)]
setorder(rd, Chromosome, BinLowEdge, BinUpEdge)

# 用各染色体数据中的最大上边界确定连续基因组坐标。
chr_layout <- rd[, .(chr_length = max(BinUpEdge, na.rm = TRUE)), by = Chromosome]
setorder(chr_layout, Chromosome)
chr_layout[, offset := shift(cumsum(chr_length), fill = 0)]
chr_layout[, center := offset + chr_length / 2]
chr_layout[, boundary := offset + chr_length]

rd <- chr_layout[rd, on = "Chromosome"]

# 使用bin中点归入固定bp窗口，并通过累计和快速计算窗口均值。
# 只输出完整的100 Kb窗口；不跨染色体，也不以缺失值补零。
make_bp_window_rd <- function(step_bp) {
  per_chr <- lapply(seq_len(nrow(chr_layout)), function(i) {
    chr_name <- as.character(chr_layout$Chromosome[i])
    chr_length <- chr_layout$chr_length[i]
    chr_offset <- chr_layout$offset[i]
    chr_rd <- rd[as.character(Chromosome) == chr_name]

    last_start <- floor((chr_length - window_bp) / step_bp) * step_bp
    starts <- seq(0, last_start, by = step_bp)
    ends <- starts + window_bp
    positions <- (chr_rd$BinLowEdge + chr_rd$BinUpEdge) / 2
    left_idx <- findInterval(starts, positions) + 1L
    right_idx <- findInterval(ends - 1e-7, positions)

    result <- data.table(
      Chromosome = factor(chr_name, levels = autosomes, ordered = TRUE),
      genome_x = chr_offset + starts + window_bp / 2
    )

    for (sample_name in sample_names) {
      values <- chr_rd[[sample_name]]
      valid <- is.finite(values)
      cumulative_sum <- c(0, cumsum(fifelse(valid, values, 0)))
      cumulative_n <- c(0L, cumsum(valid))
      window_sum <- cumulative_sum[right_idx + 1L] - cumulative_sum[left_idx]
      window_n <- cumulative_n[right_idx + 1L] - cumulative_n[left_idx]
      result[, (sample_name) := fifelse(window_n > 0L, window_sum / window_n, NA_real_)]
    }
    result
  })

  plot_rd <- rbindlist(per_chr)
  plot_rd <- melt(
    plot_rd,
    id.vars = c("Chromosome", "genome_x"),
    measure.vars = sample_names,
    variable.name = "sample",
    value.name = "read_depth",
    variable.factor = FALSE
  )
  plot_rd[, sample := factor(sample, levels = sample_names)]
  plot_rd[is.finite(read_depth)]
}

# 读取CNV并将每种样本对映射到固定轨道 ---------------------------------------
cnv <- fread(cnv_file, encoding = "UTF-8", showProgress = FALSE)
required_cnv_cols <- c(
  "pair_sampleA_name", "pair_sampleB_name", "chr", "start_pos", "end_pos"
)
missing_cnv_cols <- setdiff(required_cnv_cols, names(cnv))
if (length(missing_cnv_cols) > 0L) {
  stop("CNVresult_filter.csv缺少列: ", paste(missing_cnv_cols, collapse = ", "),
       call. = FALSE)
}

sample_rank <- setNames(seq_along(sample_names), sample_names)
canonical_pair <- function(a, b) {
  if (is.na(sample_rank[[a]]) || is.na(sample_rank[[b]])) return(NA_character_)
  ordered <- c(a, b)[order(c(sample_rank[[a]], sample_rank[[b]]))]
  paste(ordered, collapse = "-vs-")
}

cnv[, pair := mapply(canonical_pair, pair_sampleA_name, pair_sampleB_name,
                     USE.NAMES = FALSE)]
pair_levels <- c("FASR-vs-MCF7", "FASR-vs-TAMR", "MCF7-vs-TAMR")
pair_labels <- c(
  "FASR-vs-MCF7" = "FASR vs MCF7",
  "FASR-vs-TAMR" = "FASR vs TAMR",
  "MCF7-vs-TAMR" = "MCF7 vs TAMR"
)
cnv <- cnv[chr %chin% autosomes & pair %chin% pair_levels]
cnv[, Chromosome := factor(chr, levels = autosomes, ordered = TRUE)]
cnv <- chr_layout[cnv, on = "Chromosome"]
cnv[, `:=`(
  genome_start = offset + pmax(0, start_pos),
  genome_end = offset + pmin(chr_length, end_pos)
)]

cnv[, pair_display := factor(pair_labels[pair], levels = unname(pair_labels))]

# 根据保存的scaffold计算未评估区间。灰色只表示没有scaffold，不表示非CNV。 ----
unassessed <- data.table()
if (show_unassessed && file.exists(scaffold_file)) {
  scaffold <- fread(scaffold_file, header = FALSE, showProgress = FALSE)
  if (ncol(scaffold) < 2L) {
    warning("CNVscaffold_background.csv少于两列，跳过未评估区间标记。",
            call. = FALSE)
  } else {
    setnames(scaffold, names(scaffold)[1:2], c("start_idx", "end_idx"))
    valid_scaffold <- scaffold[start_idx >= 1L & end_idx >= start_idx]
    valid_scaffold[, `:=`(
      start_match = match(start_idx, rd$source_row),
      end_match = match(end_idx, rd$source_row)
    )]
    valid_scaffold <- valid_scaffold[!is.na(start_match) & !is.na(end_match)]

    # CNVscaffold中的索引对应metaDataDF的原始全局行号；使用source_row映射，
    # 避免排序或过滤后按当前行号误取坐标。
    covered <- valid_scaffold[, .(
      Chromosome = rd$Chromosome[start_match],
      end_chr = rd$Chromosome[end_match],
      start = rd$BinLowEdge[start_match],
      end = rd$BinUpEdge[end_match]
    )]
    covered <- covered[Chromosome == end_chr & !is.na(Chromosome)]
    covered[, end_chr := NULL]

    complement_intervals <- function(chr_name, chr_length_value, covered_chr) {
      if (nrow(covered_chr) == 0L) {
        return(data.table(Chromosome = chr_name, start = 0, end = chr_length_value))
      }
      setorder(covered_chr, start, end)
      covered_chr[, running_end := cummax(end)]
      covered_chr[, new_group := start > shift(running_end, fill = -Inf)]
      covered_chr[, group_id := cumsum(new_group)]
      merged <- covered_chr[, .(start = min(start), end = max(end)), by = group_id]
      starts <- c(0, merged$end)
      ends <- c(merged$start, chr_length_value)
      gaps <- data.table(start = starts, end = ends)[end > start]
      gaps[, Chromosome := chr_name]
      gaps[, .(Chromosome, start, end)]
    }

    unassessed <- rbindlist(lapply(seq_len(nrow(chr_layout)), function(i) {
      chr_name <- as.character(chr_layout$Chromosome[i])
      complement_intervals(
        chr_name,
        chr_layout$chr_length[i],
        copy(covered[as.character(Chromosome) == chr_name, .(start, end)])
      )
    }))
    unassessed[, Chromosome := factor(Chromosome, levels = autosomes, ordered = TRUE)]
    unassessed <- chr_layout[unassessed, on = "Chromosome"]
    unassessed[, `:=`(
      genome_start = offset + start,
      genome_end = offset + end,
      region_type = "Not Assessed"
    )]
  }
} else if (show_unassessed) {
  warning("未找到CNVscaffold_background.csv，图中不显示未评估区间。",
          call. = FALSE)
}

if (!show_unassessed) {
  unassessed <- data.table(
    genome_start = numeric(),
    genome_end = numeric(),
    region_type = character()
  )
}

# 绘图 -----------------------------------------------------------------------
sample_colors <- c(FASR = "#155075", MCF7 = "#5B8E3E", TAMR = "#D07C32")
pair_colors <- c(
  "FASR vs MCF7" = "#7B61A8",
  "FASR vs TAMR" = "#D64F4F",
  "MCF7 vs TAMR" = "#2A9D8F"
)
fill_colors <- pair_colors
total_length <- sum(chr_layout$chr_length)
separator_x <- head(chr_layout$boundary, -1L)

# 每个子图只保留该样本对的两条RD曲线。
pair_samples <- data.table(
  pair_display = factor(
    rep(unname(pair_labels), each = 2L),
    levels = unname(pair_labels)
  ),
  sample = factor(
    c("FASR", "MCF7", "FASR", "TAMR", "MCF7", "TAMR"),
    levels = sample_names
  )
)

draw_and_save <- function(step_bp) {
  plot_rd <- make_bp_window_rd(step_bp)
  panel_rd <- merge(plot_rd, pair_samples, by = "sample", allow.cartesian = TRUE)
  step_kb <- as.integer(step_bp / 1000)
  file_stem <- sprintf("pairwise_three_panels_window100Kb_step%dKb_no_gray", step_kb)
  output_pdf <- file.path(output_dir, paste0(file_stem, ".pdf"))
  output_png <- file.path(output_dir, paste0(file_stem, ".png"))

  p <- ggplot() +
    geom_rect(
      data = unassessed,
      aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf,
          fill = region_type),
      inherit.aes = FALSE, alpha = 0.42
    ) +
    geom_rect(
      data = cnv,
      aes(xmin = genome_start, xmax = genome_end,
          ymin = -Inf, ymax = Inf, fill = pair_display),
      inherit.aes = FALSE, alpha = 0.20
    ) +
    geom_line(
      data = panel_rd,
      aes(x = genome_x, y = read_depth, color = sample,
          group = interaction(sample, Chromosome)),
      linewidth = 0.28, alpha = 0.95, na.rm = TRUE
    ) +
    geom_vline(xintercept = separator_x, linetype = "dashed",
               linewidth = 0.3, color = "#777777") +
    facet_grid(pair_display ~ ., switch = "y") +
    scale_color_manual(values = sample_colors, name = "Read depth") +
    scale_fill_manual(values = fill_colors, name = "CNV / assessment status",
                      drop = FALSE) +
    scale_x_continuous(
      breaks = chr_layout$center,
      labels = as.character(chr_layout$Chromosome),
      limits = c(0, total_length),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = c(0, 50, 100, 150, 200, 250, 300),
      labels = label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(ylim = c(10, 300), clip = "on") +
    labs(x = NULL, y = "Mean read depth (100 Kb window)") +
    guides(
      color = guide_legend(order = 1, nrow = 1, byrow = TRUE,
                           title.position = "left", label.position = "right",
                           override.aes = list(linewidth = 1.2)),
      fill = guide_legend(order = 2, nrow = 1, byrow = TRUE,
                          title.position = "left", label.position = "right")
    ) +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(
      axis.text.x = element_text(size = 7, angle = 60, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.3),
      panel.grid = element_blank(),
      panel.spacing.y = grid::unit(0.12, "cm"),
      strip.placement = "outside",
      strip.background = element_rect(fill = "#F2F2F2", color = "#BDBDBD",
                                      linewidth = 0.25),
      strip.text.y.left = element_text(size = 8, face = "bold", angle = 90),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7.5),
      legend.key.height = grid::unit(0.34, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    )

  ggsave(output_pdf, p, width = 170, height = 110, device = cairo_pdf,
         units = "mm", family = "Arial", bg = "white")
  ggsave(output_png, p, device = grDevices::png,
         width = 170, height = 110, dpi = 600,
         units = "mm", bg = "white")

  cat("Saved step ", step_kb, " Kb outputs:\n", sep = "")
  cat("  ", output_pdf, "\n", sep = "")
  cat("  ", output_png, "\n", sep = "")
  cat("  Panel RD points: ", format(nrow(panel_rd), big.mark = ","), "\n", sep = "")
  invisible(list(plot = p, panel_rd = panel_rd))
}

results <- lapply(step_values, draw_and_save)
cat("CNV intervals: ", format(nrow(cnv), big.mark = ","), "\n", sep = "")
cat("Unassessed intervals: ", format(nrow(unassessed), big.mark = ","), "\n", sep = "")







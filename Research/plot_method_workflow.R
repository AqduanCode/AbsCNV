if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

workflow_subdir <- function(current, historical) {
  if (dir.exists(file.path(root, current))) return(current)
  if (dir.exists(file.path(root, historical))) return(historical)
  stop("Missing workflow input directory: ", current, " or ", historical)
}


samples <- c("Sample_A", "Sample_B", "Sample_C")
sample_cols <- c(Sample_A = "#3B6FB6", Sample_B = "#D97932", Sample_C = "#3B9B75")
state_cols <- c(`-1` = "#D9D9D9", `0` = "#111111", `1` = "#4C78A8",
                `2` = "#A9A9A9", `3` = "#F2A541", `4` = "#C44E52")
pair_levels <- c("A_vs_B", "A_vs_C", "B_vs_C")
pair_labels <- c(A_vs_B = "Sample_A vs Sample_B",
                 A_vs_C = "Sample_A vs Sample_C",
                 B_vs_C = "Sample_B vs Sample_C")

# All visible text is at least 6 pt. The y-axis title receives extra right margin.
paper_theme <- theme_bw(base_size = 7.2) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.18, colour = "#E7E7E7"),
    axis.text = element_text(size = 6.2, colour = "#333333"),
    axis.title = element_text(size = 7.2, colour = "#222222"),
    axis.title.y = element_text(margin = margin(r = 9)),
    axis.title.x = element_text(margin = margin(t = 3)),
    strip.text = element_text(size = 6.4, colour = "#222222"),
    strip.background = element_rect(fill = "#ECECEC", colour = "#555555", linewidth = 0.25),
    legend.text = element_text(size = 6.0),
    legend.title = element_text(size = 6.4),
    legend.key.height = unit(2.8, "mm"),
    legend.key.width = unit(4.5, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    plot.margin = margin(4, 5, 4, 7)
  )

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 8.5, colour = "#111111"),
  plot.tag.position = c(0.002, 0.995)
)

norm <- fread(file.path(root, workflow_subdir("01_Fig1_RD_normalization_smoothing", "01_Fig1b_RD_normalization_smoothing"),
                        "normalized_and_smooth10_RD.csv.gz"))
params <- fread(file.path(root, workflow_subdir("02_Fig1_sample_specific_GMM", "02_Fig1c_sample_specific_GMM"),
                          "gmm_parameters_from_fitting_gmm.csv"))
bp <- fread(file.path(root, workflow_subdir("03_Fig1_breakpoints_shared_scaffold", "03_Fig1d_breakpoints_shared_scaffold"),
                      "breakpoint_scores_from_AbsCNV.csv.gz"))
scaffold <- fread(file.path(root, workflow_subdir("03_Fig1_breakpoints_shared_scaffold", "03_Fig1d_breakpoints_shared_scaffold"),
                            "shared_scaffold_from_genome_segment.csv"))
merged_runs <- fread(file.path(root, workflow_subdir("03_Fig1_breakpoints_shared_scaffold", "03_Fig1d_breakpoints_shared_scaffold"),
                               "merged_breakpoint_runs_from_genome_segment_logic.csv"))
calls <- fread(file.path(root, workflow_subdir("04_Fig1_posterior_state_calls", "04_Fig1e_posterior_state_calls"),
                         "shared_segments_and_state_calls.csv"))
pairs <- fread(file.path(root, workflow_subdir("05_Fig1_pairwise_refinement", "05_Fig1f_pairwise_refinement"),
                         "pairwise_segment_calls_from_AbsCNV.csv"))

# b: normalized read depth and 10-bin smoothing.
norm_long <- rbindlist(lapply(samples, function(s) {
  data.table(position_mb = norm$binLowEdge / 1e6, sample = s,
             normalized_rd = norm[[paste0(s, "_normalized_rd")]],
             smooth10_rd = norm[[paste0(s, "_smooth10_rd")]])
}))
p_b <- ggplot(norm_long, aes(position_mb)) +
  geom_point(aes(y = normalized_rd), size = 0.10, alpha = 0.14, colour = "#8F8F8F") +
  geom_line(aes(y = smooth10_rd, colour = sample), linewidth = 0.25, na.rm = TRUE) +
  geom_hline(yintercept = 100, linetype = 3, linewidth = 0.25, colour = "#555555") +
  facet_wrap(~sample, ncol = 1) +
  scale_colour_manual(values = sample_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = c(50, 100, 150, 200)) +
  coord_cartesian(ylim = c(20, 225), clip = "off") +
  labs(tag = "b", x = "Position on chrSim (Mb)", y = "Normalized RD") +
  paper_theme + tag_theme +
  theme(panel.spacing.y = unit(0.7, "mm"),
        strip.text = element_text(margin = margin(0.7, 1, 0.7, 1)),
        plot.margin = margin(4, 5, 2, 8))

# c: sample-specific GMM.
density_data <- rbindlist(lapply(samples, function(s) {
  x <- norm[[paste0(s, "_smooth10_rd")]]
  data.table(sample = s, value = x[is.finite(x)])
}))
curve_data <- rbindlist(lapply(samples, function(s) {
  ps <- params[sample == s]
  x <- seq(20, 225, length.out = 900)
  rbindlist(lapply(seq_len(nrow(ps)), function(i) {
    data.table(sample = s, state = factor(ps$state[i]), x = x,
               density = ps$weight[i] * dnorm(x, ps$mean[i], ps$std[i]))
  }))
}))
p_c <- ggplot(density_data, aes(value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 fill = "#D9D9D9", colour = NA) +
  geom_line(data = curve_data, aes(x = x, y = density, colour = state), linewidth = 0.48) +
  facet_wrap(~sample, nrow = 1) +
  scale_colour_manual(values = state_cols[c("1", "2", "3", "4")]) +
  labs(tag = "c", x = "10-bin mean normalized RD", y = "Density",
       colour = "CN-like state") +
  paper_theme + tag_theme +
  theme(panel.grid = element_blank(), legend.position = "top",
        legend.direction = "horizontal", legend.margin = margin(0, 0, -2, 0),
        plot.margin = margin(3, 5, 2, 8)) +
  guides(colour = guide_legend(nrow = 1, title.position = "left",
                               override.aes = list(linewidth = 0.8)))

# d: RD tracks, candidate breakpoints, merged runs and shared scaffold.
rd_tracks <- rbindlist(lapply(samples, function(s) {
  data.table(position_mb = norm$binLowEdge / 1e6, sample = s,
             smooth10_rd = norm[[paste0(s, "_smooth10_rd")]])
}))
sample_breaks <- rbindlist(lapply(samples, function(s) {
  data.table(position_mb = bp$binLowEdge[bp[[paste0(s, "_candidate")]] == 1] / 1e6,
             sample = s)
}))
sample_breaks[, track_y := c(Sample_A = 3, Sample_B = 2, Sample_C = 1)[sample]]
boundary_lines <- sort(unique(c(scaffold$start_bp, scaffold$end_bp))) / 1e6
scaffold_plot <- copy(scaffold)
scaffold_plot[, `:=`(
  xmin = pmin(start_bp / 1e6 + 0.01, (start_bp + end_bp) / 2e6),
  xmax = pmax(end_bp / 1e6 - 0.01, (start_bp + end_bp) / 2e6),
  scaffold_fill = ifelse((segment_id - 1L) %% 2L == 0L, "dark", "light")
)]
p_d_top <- ggplot() +
  geom_line(data = rd_tracks, aes(position_mb, smooth10_rd, colour = sample),
            linewidth = 0.22, na.rm = TRUE) +
  geom_vline(data = sample_breaks, aes(xintercept = position_mb, colour = sample),
             linewidth = 0.16, alpha = 0.48) +
  facet_grid(rows = vars(sample), scales = "free_y") +
  scale_colour_manual(values = sample_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = c(50, 100, 150, 200)) +
  labs(y = "10-bin mean RD", x = NULL) +
  paper_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing.y = unit(0.6, "mm"),
        strip.text.y = element_text(angle = 0),
        plot.margin = margin(4, 5, 0, 8))
p_d_bottom <- ggplot() +
  geom_rect(data = scaffold_plot,
            aes(xmin = xmin, xmax = xmax, ymin = -0.92, ymax = -0.34,
                fill = scaffold_fill), colour = "#7A3E9D", linewidth = 0.22) +
  geom_rect(data = merged_runs,
            aes(xmin = start_bp / 1e6, xmax = end_bp / 1e6,
                ymin = 0.12, ymax = 0.62), fill = "#C9B2DA", colour = NA, alpha = 0.48) +
  geom_segment(data = sample_breaks,
               aes(x = position_mb, xend = position_mb,
                   y = track_y - 0.30, yend = track_y + 0.30, colour = sample),
               linewidth = 0.24, alpha = 0.78) +
  geom_segment(data = data.table(position_mb = boundary_lines),
               aes(x = position_mb, xend = position_mb, y = -0.98, yend = -0.28),
               linewidth = 0.30, colour = "#7A3E9D") +
  scale_colour_manual(values = sample_cols, guide = "none") +
  scale_fill_manual(values = c(dark = "#666666", light = "#B8B8B8"), guide = "none") +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = c(3, 2, 1, 0.37, -0.63),
                     labels = c("Sample_A", "Sample_B", "Sample_C",
                                "Merged run", "Shared scaffold"),
                     limits = c(-1.02, 3.45), position = "right") +
  labs(x = "Position on chrSim (Mb)", y = NULL) +
  paper_theme + theme(plot.margin = margin(0, 5, 2, 8))

draw_d <- function(vp = NULL, tag = TRUE) {
  if (!is.null(vp)) pushViewport(vp)
  aligned_d <- cowplot::align_plots(p_d_top, p_d_bottom,
                                    align = "v", axis = "lr")
  lay <- grid.layout(2, 1, heights = unit(c(27, 28), "null"))
  pushViewport(viewport(layout = lay))
  pushViewport(viewport(layout.pos.row = 1))
  grid.draw(aligned_d[[1]])
  upViewport()
  pushViewport(viewport(layout.pos.row = 2))
  grid.draw(aligned_d[[2]])
  upViewport()
  upViewport()
  if (tag) grid.text("d", x = unit(1.3, "mm"), y = unit(1, "npc") - unit(1.2, "mm"),
                     just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 8.5))
  if (!is.null(vp)) upViewport()
}

# e: posterior state calls.
call_long <- rbindlist(lapply(samples, function(s) {
  data.table(segment_id = calls$segment_id, sample = s,
             called_state = calls[[paste0(s, "_called_state")]],
             max_posterior = calls[[paste0(s, "_max_posterior")]])
}))
call_long[, label := ifelse(called_state == -1, "NC",
                            paste0(called_state, "\nP=", sprintf("%.2f", max_posterior)))]
p_e <- ggplot(call_long,
              aes(factor(segment_id), factor(sample, levels = rev(samples)),
                  fill = factor(called_state))) +
  geom_tile(colour = "white", linewidth = 0.18) +
  geom_text(aes(label = label), size = 2.2, lineheight = 0.88) +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  labs(tag = "e", x = "Shared segment ID", y = NULL, fill = "Called state") +
  paper_theme + tag_theme +
  theme(panel.grid = element_blank(), legend.position = "top",
        legend.margin = margin(0, 0, -2, 0),
        plot.margin = margin(3, 5, 2, 8)) +
  guides(fill = guide_legend(nrow = 1, title.position = "left"))

# f: pairwise refinement.
pairs[, pair := factor(pair, levels = pair_levels, labels = unname(pair_labels[pair_levels]))]
pairs[, `:=`(xmin = start_bp / 1e6, xmax = end_bp / 1e6)]
pairs[, direction := fifelse(
  refined_positive == 0, "No pairwise CNV",
  fifelse(sample_1_state < 0 | sample_2_state < 0, "CNV with no-call",
          fifelse(sample_1_state > sample_2_state, "First sample higher CN",
                  fifelse(sample_2_state > sample_1_state, "Second sample higher CN",
                          "Different RD; same called CN"))))]
pairs[, cn_label := ifelse(refined_positive == 1,
                           paste0("CN", sample_1_state, " vs CN", sample_2_state), "")]
pair_members <- list(
  "Sample_A vs Sample_B" = c("Sample_A", "Sample_B"),
  "Sample_A vs Sample_C" = c("Sample_A", "Sample_C"),
  "Sample_B vs Sample_C" = c("Sample_B", "Sample_C"))
pair_rd <- rbindlist(lapply(names(pair_members), function(pair_name) {
  rbindlist(lapply(pair_members[[pair_name]], function(s) {
    data.table(pair = factor(pair_name, levels = unname(pair_labels[pair_levels])),
               sample = s, position_mb = norm$binLowEdge / 1e6,
               smooth10_rd = norm[[paste0(s, "_smooth10_rd")]])
  }))
}))
pair_cols <- c("No pairwise CNV" = "#E3E3E3",
               "First sample higher CN" = "#3B6FB6",
               "Second sample higher CN" = "#D97932",
               "Different RD; same called CN" = "#7A3E9D",
               "CNV with no-call" = "#C44E52")
p_f <- ggplot(pairs) +
  geom_rect(data = pairs[refined_positive == 1],
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = direction),
            inherit.aes = FALSE, alpha = 0.18, colour = NA) +
  geom_line(data = pair_rd, aes(position_mb, smooth10_rd, colour = sample),
            linewidth = 0.22, na.rm = TRUE) +
  geom_text(data = pairs[refined_positive == 1 & (xmax - xmin) >= 0.30],
            aes(x = (xmin + xmax) / 2, y = 218, label = cn_label),
            inherit.aes = FALSE, size = 2.2, check_overlap = TRUE) +
  facet_grid(rows = vars(pair)) +
  scale_fill_manual(values = pair_cols, drop = TRUE) +
  scale_colour_manual(values = sample_cols, name = "10-bin mean RD") +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = c(50, 100, 150, 200)) +
  coord_cartesian(ylim = c(20, 225), clip = "off") +
  labs(tag = "f", x = "Position on chrSim (Mb)", y = "10-bin mean RD",
       fill = "Refined pairwise CNV") +
  paper_theme + tag_theme +
  theme(panel.spacing.y = unit(0.7, "mm"), legend.position = "top",
        legend.box = "vertical", legend.margin = margin(0, 0, -2, 0),
        strip.text.y = element_text(angle = 0),
        plot.margin = margin(3, 5, 2, 8)) +
  guides(colour = guide_legend(nrow = 1, title.position = "left",
                               override.aes = list(linewidth = 0.65)),
         fill = guide_legend(nrow = 1, title.position = "left"))

mm_to_in <- function(x) x / 25.4
save_one <- function(plot, name, height_mm) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = mm_to_in(170), height = mm_to_in(height_mm),
         device = cairo_pdf, bg = "white", limitsize = FALSE)
}
save_one(p_b, "Fig1b_redrawn_170mm", 36)
save_one(p_c, "Fig1c_redrawn_170mm", 30)
save_one(p_e, "Fig1e_redrawn_170mm", 28)
save_one(p_f, "Fig1f_redrawn_170mm", 46)
cairo_pdf(file.path(out_dir, "Fig1d_redrawn_170mm.pdf"),
          width = mm_to_in(170), height = mm_to_in(55), bg = "white")
grid.newpage(); draw_d(tag = FALSE); dev.off()

# Combined b-f figure: 170 x 195 mm.
draw_combined <- function() {
  grid.newpage()
  lay <- grid.layout(5, 1, heights = unit(c(36, 30, 55, 28, 46), "mm"))
  pushViewport(viewport(layout = lay))
  print(p_b, vp = viewport(layout.pos.row = 1), newpage = FALSE)
  print(p_c, vp = viewport(layout.pos.row = 2), newpage = FALSE)
  draw_d(viewport(layout.pos.row = 3))
  print(p_e, vp = viewport(layout.pos.row = 4), newpage = FALSE)
  print(p_f, vp = viewport(layout.pos.row = 5), newpage = FALSE)
  upViewport()
}

cairo_pdf(file.path(out_dir, "Fig1_b-f_combined_170x195mm.pdf"),
          width = mm_to_in(170), height = mm_to_in(195), bg = "white")
draw_combined(); dev.off()

png(file.path(out_dir, "Fig1_b-f_combined_170x195mm.png"),
    width = round(mm_to_in(170) * 300), height = round(mm_to_in(195) * 300),
    res = 300, bg = "white")
draw_combined(); dev.off()

svg(file.path(out_dir, "Fig1_b-f_combined_170x195mm.svg"),
    width = mm_to_in(170), height = mm_to_in(195), onefile = TRUE, bg = "white")
draw_combined(); dev.off()

cat("Redrawn panels and combined figure written to:", out_dir, "\n")

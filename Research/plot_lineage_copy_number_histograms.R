# Reproduce the moving-average CN-peak plots for the background-updated results.
# The 10 sample panels are combined into one 5-column x 2-row PDF.

options(warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese"))
  windowsFonts(Arial = windowsFont("Arial"))
}

project_env <- Sys.getenv("CNV_PLOT_PROJECT_DIR", unset = "")
output_env <- Sys.getenv("CNV_PLOT_OUTPUT_DIR", unset = "")
if (nzchar(project_env) && nzchar(output_env)) {
  project_dir <- project_env
  output_dir <- output_env
} else {
  script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_dir <- if (length(script_argument)) {
    dirname(normalizePath(sub("^--file=", "", script_argument[1]),
                          winslash = "/", mustWork = TRUE))
  } else {
    normalizePath(".", winslash = "/", mustWork = TRUE)
  }
  project_dir <- dirname(script_dir)
  output_dir <- script_dir
}
setwd(project_dir)

metadata_file <- "output/metaDataDF.csv"
parameter_file <- file.path(output_dir, "cn_peak_parameters_background_updated.csv")
output_pdf <- file.path(output_dir, "all_samples_moving_average_histograms_5x2.pdf")
summary_file <- file.path(output_dir, "moving_average_histogram_summary.csv")

window_size <- 50L
bin_width <- 2
x_limits <- c(0, 250)
figure_width_mm <- 175
figure_height_mm <- 92

sample_map <- c(
  G6 = "U1",
  G6A7first = "U1.1",
  G6B11first = "U1.2",
  G6C7first = "U1.3",
  G6A7E1 = "U1.1.1",
  G6A7E11 = "U1.1.2",
  G6B11C1 = "U1.2.1",
  G6B11C6 = "U1.2.2",
  G6C7E1 = "U1.3.1",
  G6C7E9 = "U1.3.2"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(metadata_file)) stop("Missing metadata file: ", metadata_file)
if (!file.exists(parameter_file)) stop("Missing peak parameter file: ", parameter_file)

metadata_columns <- names(fread(metadata_file, nrows = 0L))
raw_samples <- metadata_columns[-seq_len(3L)]
if (!identical(raw_samples, names(sample_map))) {
  stop("Unexpected sample columns/order in metaDataDF.csv: ",
       paste(raw_samples, collapse = ", "))
}

peak_parameters <- fread(parameter_file)
required_peak_columns <- c("raw_sample_name", "sample_name", "copy_number", "mean")
if (!all(required_peak_columns %in% names(peak_parameters))) {
  stop("Peak parameter table must contain: ", paste(required_peak_columns, collapse = ", "))
}

message("Reading ", length(raw_samples), " samples from: ", metadata_file)
depth_data <- fread(metadata_file, select = raw_samples, showProgress = TRUE)

histogram_rows <- vector("list", length(raw_samples))
peak_rows <- vector("list", length(raw_samples))
summary_rows <- vector("list", length(raw_samples))

for (index in seq_along(raw_samples)) {
  raw_name <- raw_samples[index]
  display_name <- unname(sample_map[[raw_name]])
  message("Processing [", index, "/", length(raw_samples), "]: ", raw_name,
          " -> ", display_name)

  depth <- as.numeric(depth_data[[raw_name]])
  depth <- depth[is.finite(depth)]
  if (length(depth) < window_size) stop("Too few observations for ", raw_name)

  moving_average <- frollmean(depth, n = window_size, align = "left", na.rm = FALSE)
  moving_average <- moving_average[is.finite(moving_average)]
  breaks <- seq(x_limits[1], x_limits[2], by = bin_width)
  histogram <- hist(moving_average, breaks = breaks, plot = FALSE,
                    right = TRUE, include.lowest = TRUE)
  histogram_rows[[index]] <- data.table(
    sample_name = display_name,
    x = histogram$mids,
    count = histogram$counts
  )

  this_peaks <- copy(peak_parameters[raw_sample_name == raw_name][order(copy_number)])
  if (!nrow(this_peaks)) stop("No peak parameters found for ", raw_name)
  this_peaks[, peak_count := vapply(mean, function(mu) {
    local_bins <- which(abs(histogram$mids - mu) <= 3 * bin_width)
    if (!length(local_bins)) return(0)
    max(histogram$counts[local_bins])
  }, numeric(1))]
  # Two-line labels preserve the old CN-to-peak-mean information in narrow facets.
  this_peaks[, label := sprintf("CN%d\n%.1f", copy_number, mean)]
  this_peaks[, sample_name := display_name]
  peak_rows[[index]] <- this_peaks

  summary_rows[[index]] <- data.table(
    raw_sample_name = raw_name,
    sample_name = display_name,
    input_bins = length(depth),
    moving_average_windows = length(moving_average),
    moving_average_mean = mean(moving_average),
    moving_average_sd = sd(moving_average)
  )
}

histogram_data <- rbindlist(histogram_rows)
plot_peaks <- rbindlist(peak_rows)
summary_data <- rbindlist(summary_rows)

display_order <- unname(sample_map)
histogram_data[, sample_name := factor(sample_name, levels = display_order)]
plot_peaks[, sample_name := factor(sample_name, levels = display_order)]

global_y_max <- max(histogram_data$count, na.rm = TRUE)
sample_labels <- data.table(
  sample_name = factor(display_order, levels = display_order),
  x = 7,
  y = global_y_max * 1.155,
  label = display_order
)
plot_peaks[, label_y := pmax(peak_count * 1.025, global_y_max * 0.08)]

plot_object <- ggplot(histogram_data, aes(x = x, y = count)) +
  geom_col(width = bin_width, fill = "#9ECAE1", color = "#6BAED6",
           linewidth = 0.10, alpha = 0.85) +
  geom_segment(
    data = plot_peaks,
    aes(x = mean, xend = mean, y = 0, yend = peak_count),
    inherit.aes = FALSE, color = "#CB181D", linewidth = 0.38
  ) +
  geom_text(
    data = plot_peaks,
    aes(x = mean, y = label_y, label = label),
    inherit.aes = FALSE, family = "Arial", size = 5.8 / ggplot2::.pt,
    color = "#A50F15", fontface = "bold", lineheight = 0.82, vjust = 0
  ) +
  geom_text(
    data = sample_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, family = "Arial", size = 7.5 / ggplot2::.pt,
    fontface = "bold", hjust = 0, vjust = 1
  ) +
  facet_wrap(~sample_name, ncol = 5, nrow = 2, drop = FALSE) +
  scale_x_continuous(limits = x_limits, breaks = c(0, 100, 200),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, global_y_max * 1.18),
                     breaks = scales::breaks_pretty(n = 3),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "50-bin moving-average read depth", y = "Frequency") +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 6.2, color = "black"),
    axis.title = element_text(size = 7.2, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.30, color = "black"),
    panel.spacing.x = unit(2.2, "mm"),
    panel.spacing.y = unit(3.5, "mm"),
    plot.margin = margin(2.5, 2.5, 2.5, 2.5, unit = "mm")
  )

ggsave(output_pdf, plot_object, width = figure_width_mm, height = figure_height_mm,
       units = "mm", device = cairo_pdf, family = "Arial", bg = "white")
fwrite(summary_data, summary_file)

message("Completed combined figure: ", output_pdf)
message("Figure dimensions: ", figure_width_mm, " x ", figure_height_mm, " mm")


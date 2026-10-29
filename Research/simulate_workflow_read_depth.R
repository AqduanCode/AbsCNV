if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], winslash = "/", mustWork = FALSE) else normalizePath("..", winslash = "/", mustWork = FALSE)
noise_profile <- if (length(args) >= 2) args[2] else "revised"
out_dir <- file.path(root, "00_simulation_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260725L
set.seed(seed)
bin_size <- 1000L
n_bins <- 20000L
samples <- c("Sample_A", "Sample_B", "Sample_C")
pos <- 0:(n_bins - 1L)

set_state <- function(x, start_mb, end_mb, state) {
  start_idx <- as.integer(start_mb * 1000) + 1L
  end_idx <- as.integer(end_mb * 1000)
  x[start_idx:end_idx] <- state
  x
}

truth <- setNames(lapply(samples, function(x) rep(2, n_bins)), samples)

# Shared state anchors provide observations for GMM states 1, 3 and 4.
for (s in samples) {
  truth[[s]] <- set_state(truth[[s]], 0.50, 1.00, 1)
  truth[[s]] <- set_state(truth[[s]], 1.50, 2.00, 3)
  truth[[s]] <- set_state(truth[[s]], 2.50, 3.00, 4)
}

# Displayed CNV-like events. Slightly shifted boundaries test shared scaffolding.
truth[["Sample_A"]] <- set_state(truth[["Sample_A"]], 4.00, 6.00, 1)
truth[["Sample_B"]] <- set_state(truth[["Sample_B"]], 4.05, 6.00, 1)
truth[["Sample_C"]] <- set_state(truth[["Sample_C"]], 4.00, 6.10, 1)

truth[["Sample_A"]] <- set_state(truth[["Sample_A"]], 8.00, 10.00, 3)
truth[["Sample_C"]] <- set_state(truth[["Sample_C"]], 8.00, 10.00, 3)

truth[["Sample_B"]] <- set_state(truth[["Sample_B"]], 12.00, 14.00, 4)
truth[["Sample_C"]] <- set_state(truth[["Sample_C"]], 12.05, 13.95, 4)

# Ambiguous effective state in A; C supplies clear shared boundaries.
truth[["Sample_A"]] <- set_state(truth[["Sample_A"]], 15.50, 16.20, 2.40)
truth[["Sample_C"]] <- set_state(truth[["Sample_C"]], 15.50, 16.20, 3)

# C-specific loss.
truth[["Sample_C"]] <- set_state(truth[["Sample_C"]], 17.00, 18.50, 1)

baseline <- c(Sample_A = 78, Sample_B = 104, Sample_C = 132)
if (noise_profile == "original") {
  dispersion <- c(Sample_A = 18, Sample_B = 24, Sample_C = 30)
  wave_slow <- 0.035
  wave_fast <- 0.018
} else if (noise_profile == "revised") {
  # Keep different sample depths, but use the same biological-dispersion
  # parameter so the low-depth sample is not uniquely penalized.
  dispersion <- c(Sample_A = 45, Sample_B = 45, Sample_C = 45)
  wave_slow <- 0.020
  wave_fast <- 0.010
} else if (noise_profile == "low_noise_test") {
  dispersion <- c(Sample_A = 80, Sample_B = 80, Sample_C = 80)
  wave_slow <- 0.012
  wave_fast <- 0.006
} else {
  stop("Unknown noise profile: ", noise_profile)
}
gc_wave <- 1 + wave_slow * sin(2 * pi * pos / 1700) +
  wave_fast * cos(2 * pi * pos / 430)

raw <- data.table(
  chr = rep("chrSim", n_bins),
  binLowEdge = pos * bin_size,
  binUpEdge = (pos + 1L) * bin_size
)
truth_bins <- copy(raw)

for (s in samples) {
  mu <- baseline[[s]] * truth[[s]] / 2 * gc_wave
  size <- dispersion[[s]]
  counts <- rnbinom(n_bins, size = size, mu = mu)
  raw[[s]] <- counts
  truth_bins[[paste0(s, "_truth_effective_state")]] <- truth[[s]]
}

events <- rbindlist(list(
  data.table(event_id = "ANCHOR_CN1", sample = samples, start_bp = 0.50e6, end_bp = 1.00e6, truth_state = 1, role = "GMM_anchor"),
  data.table(event_id = "ANCHOR_CN3", sample = samples, start_bp = 1.50e6, end_bp = 2.00e6, truth_state = 3, role = "GMM_anchor"),
  data.table(event_id = "ANCHOR_CN4", sample = samples, start_bp = 2.50e6, end_bp = 3.00e6, truth_state = 4, role = "GMM_anchor"),
  data.table(event_id = "LOSS_SHARED", sample = samples, start_bp = c(4.00, 4.05, 4.00) * 1e6, end_bp = c(6.00, 6.00, 6.10) * 1e6, truth_state = 1, role = "display_event"),
  data.table(event_id = "GAIN_AC", sample = c("Sample_A", "Sample_C"), start_bp = 8.00e6, end_bp = 10.00e6, truth_state = 3, role = "display_event"),
  data.table(event_id = "GAIN_BC", sample = c("Sample_B", "Sample_C"), start_bp = c(12.00, 12.05) * 1e6, end_bp = c(14.00, 13.95) * 1e6, truth_state = 4, role = "display_event"),
  data.table(event_id = "AMBIGUOUS_A", sample = "Sample_A", start_bp = 15.50e6, end_bp = 16.20e6, truth_state = 2.40, role = "ambiguous"),
  data.table(event_id = "GAIN_C_BOUNDARY", sample = "Sample_C", start_bp = 15.50e6, end_bp = 16.20e6, truth_state = 3, role = "display_event"),
  data.table(event_id = "LOSS_C_PRIVATE", sample = "Sample_C", start_bp = 17.00e6, end_bp = 18.50e6, truth_state = 1, role = "display_event")
), fill = TRUE)
events[, `:=`(chr = "chrSim", length_bp = end_bp - start_bp)]
setcolorder(events, c("event_id", "sample", "chr", "start_bp", "end_bp", "length_bp", "truth_state", "role"))

fwrite(raw, file.path(out_dir, "simulated_raw_RD_20Mb_3samples.csv"))
fwrite(truth_bins, file.path(out_dir, "simulation_truth_bin_level.csv.gz"))
fwrite(events, file.path(out_dir, "simulation_truth_events.csv"))

summary <- list(
  seed = seed,
  chromosome = "chrSim",
  chromosome_length_bp = n_bins * bin_size,
  bin_size = bin_size,
  n_bins = n_bins,
  samples = samples,
  baseline_depth = as.list(baseline),
  negative_binomial_size = as.list(dispersion),
  wave_slow_amplitude = wave_slow,
  wave_fast_amplitude = wave_fast,
  noise_profile = noise_profile,
  data_role = "controlled workflow validation; not a comparative benchmark"
)
write_json(summary, file.path(out_dir, "simulation_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

cat("Simulation completed:", n_bins, "bins x", length(samples), "samples\n")

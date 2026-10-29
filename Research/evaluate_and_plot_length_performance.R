if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))

argv <- commandArgs(TRUE)
if(length(argv)!=2) stop("Usage: Rscript SCRIPT LENGTH_DATASET NEW_RESULT_DIR")
project_path <- argv[1]
output_dir <- file.path(argv[2], "evaluation")
figure_dir <- file.path(argv[2], "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

labels <- c("20kb", "50kb", "100kb", "200kb", "500kb", "1000kb")
length_bp_values <- c(20000, 50000, 100000, 200000, 500000, 1000000)
threshold <- 0.80

safe_ratio <- function(numerator, denominator) {
  if (denominator == 0) return(NA_real_)
  numerator / denominator
}

f1_value <- function(precision, recall) {
  if (is.na(precision) || is.na(recall)) return(NA_real_)
  if (precision + recall == 0) return(0)
  2 * precision * recall / (precision + recall)
}

empty_matches <- function() {
  data.frame(
    prediction_index = integer(), truth_index = integer(),
    prediction_id = character(), truth_id = character(),
    overlap_bp = numeric(), prediction_overlap_fraction = numeric(),
    truth_overlap_fraction = numeric(), iou = numeric(), score = numeric()
  )
}

greedy_match <- function(prediction, truth) {
  candidates <- list()
  candidate_index <- 0L
  if (nrow(prediction) > 0 && nrow(truth) > 0) {
    for (pi in seq_len(nrow(prediction))) {
      possible <- which(
        truth$chrom == prediction$chrom[pi] &
          truth$direction == prediction$direction[pi]
      )
      for (ti in possible) {
        overlap_bp <- max(
          0,
          min(prediction$end[pi], truth$end[ti]) -
            max(prediction$start[pi], truth$start[ti])
        )
        if (overlap_bp <= 0) next
        pred_fraction <- overlap_bp / prediction$length_bp[pi]
        truth_fraction <- overlap_bp / truth$length_bp[ti]
        if (!(pred_fraction >= threshold || truth_fraction >= threshold)) next
        union_bp <- prediction$length_bp[pi] + truth$length_bp[ti] - overlap_bp
        candidate_index <- candidate_index + 1L
        candidates[[candidate_index]] <- data.frame(
          prediction_index = pi,
          truth_index = ti,
          prediction_id = prediction$prediction_id[pi],
          truth_id = truth$truth_id[ti],
          overlap_bp = overlap_bp,
          prediction_overlap_fraction = pred_fraction,
          truth_overlap_fraction = truth_fraction,
          iou = overlap_bp / union_bp,
          score = min(pred_fraction, truth_fraction)
        )
      }
    }
  }
  if (length(candidates) == 0) return(empty_matches())
  candidates <- do.call(rbind, candidates)
  candidates <- candidates[order(
    -candidates$score, -candidates$iou, -candidates$overlap_bp,
    candidates$prediction_id, candidates$truth_id
  ), , drop = FALSE]
  used_prediction <- rep(FALSE, nrow(prediction))
  used_truth <- rep(FALSE, nrow(truth))
  selected <- list()
  selected_index <- 0L
  for (i in seq_len(nrow(candidates))) {
    pi <- candidates$prediction_index[i]
    ti <- candidates$truth_index[i]
    if (used_prediction[pi] || used_truth[ti]) next
    used_prediction[pi] <- TRUE
    used_truth[ti] <- TRUE
    selected_index <- selected_index + 1L
    selected[[selected_index]] <- candidates[i, , drop = FALSE]
  }
  if (length(selected) == 0) return(empty_matches())
  do.call(rbind, selected)
}

event_rows <- list()
base_rows <- list()
all_matches <- list()
all_unmatched_predictions <- list()
all_unmatched_truth <- list()

for (i in seq_along(labels)) {
  label <- labels[i]
  nominal_length <- length_bp_values[i]
  prediction_path <- file.path(project_path, "output", label, "CNVresult_filter.csv")
  truth_path <- file.path(
    project_path, "benchmark",
    paste0(label, "_Sample_B_vs_Sample_A.benchmark.tsv")
  )
  if (!file.exists(prediction_path)) stop("Missing prediction: ", prediction_path)
  if (!file.exists(truth_path)) stop("Missing benchmark: ", truth_path)

  prediction <- read.csv(prediction_path, check.names = FALSE)
  truth <- read.delim(truth_path, check.names = FALSE)
  required_prediction <- c("chr", "start_pos", "end_pos", "copy_number_A", "copy_number_B")
  required_truth <- c("event_id", "chrom", "start", "end", "length_bp", "delta_B_minus_A")
  # The method returns an empty two-column table when there are no positive
  # segments. Treat that as zero predictions, not as a malformed result.
  if (nrow(prediction) == 0 && length(setdiff(required_prediction, names(prediction))) > 0) {
    prediction <- data.frame(
      chr = character(), start_pos = numeric(), end_pos = numeric(),
      copy_number_A = numeric(), copy_number_B = numeric()
    )
  } else if (length(setdiff(required_prediction, names(prediction))) > 0) {
    stop(label, " prediction columns missing: ", paste(setdiff(required_prediction, names(prediction)), collapse = ", "))
  }
  if (length(setdiff(required_truth, names(truth))) > 0) {
    stop(label, " truth columns missing: ", paste(setdiff(required_truth, names(truth)), collapse = ", "))
  }

  prediction$prediction_id <- sprintf("%s_MyMethod_%03d", label, seq_len(nrow(prediction)))
  prediction$chrom <- as.character(prediction$chr)
  prediction$start <- as.numeric(prediction$start_pos)
  prediction$end <- as.numeric(prediction$end_pos)
  prediction$length_bp <- prediction$end - prediction$start
  prediction$delta <- as.numeric(prediction$copy_number_B) - as.numeric(prediction$copy_number_A)
  prediction$direction <- ifelse(prediction$delta > 0, "dup", ifelse(prediction$delta < 0, "del", "neutral"))
  prediction <- prediction[prediction$direction != "neutral", , drop = FALSE]

  truth$truth_id <- as.character(truth$event_id)
  truth$chrom <- as.character(truth$chrom)
  truth$start <- as.numeric(truth$start)
  truth$end <- as.numeric(truth$end)
  truth$length_bp <- as.numeric(truth$length_bp)
  truth$delta <- as.numeric(truth$delta_B_minus_A)
  truth$direction <- ifelse(truth$delta > 0, "dup", ifelse(truth$delta < 0, "del", "neutral"))
  truth <- truth[truth$direction != "neutral", , drop = FALSE]

  if (nrow(truth) != 20) stop(label, " benchmark does not contain 20 non-neutral events")
  if (any(truth$length_bp != nominal_length)) stop(label, " benchmark contains an unexpected length")
  if (any(prediction$length_bp <= 0) || any(truth$length_bp <= 0)) stop(label, " contains a non-positive interval")

  matches <- greedy_match(prediction, truth)
  matched_prediction <- rep(FALSE, nrow(prediction))
  matched_truth <- rep(FALSE, nrow(truth))
  if (nrow(matches) > 0) {
    matched_prediction[matches$prediction_index] <- TRUE
    matched_truth[matches$truth_index] <- TRUE
    matches$length_group <- label
    matches$chrom <- prediction$chrom[matches$prediction_index]
    matches$direction <- prediction$direction[matches$prediction_index]
    matches$prediction_start <- prediction$start[matches$prediction_index]
    matches$prediction_end <- prediction$end[matches$prediction_index]
    matches$truth_start <- truth$start[matches$truth_index]
    matches$truth_end <- truth$end[matches$truth_index]
    all_matches[[label]] <- matches
  }

  tp <- nrow(matches)
  fp <- nrow(prediction) - tp
  fn <- nrow(truth) - tp
  event_precision <- safe_ratio(tp, tp + fp)
  event_recall <- safe_ratio(tp, tp + fn)
  event_rows[[label]] <- data.frame(
    length_group = label, length_bp = nominal_length,
    predictions = nrow(prediction), truth_events = nrow(truth),
    TP = tp, FP = fp, FN = fn,
    precision = event_precision, recall = event_recall,
    F1 = f1_value(event_precision, event_recall)
  )

  predicted_bp <- sum(prediction$length_bp)
  truth_bp <- sum(truth$length_bp)
  tp_bp <- if (nrow(matches) == 0) 0 else sum(matches$overlap_bp)
  fp_bp <- predicted_bp - tp_bp
  fn_bp <- truth_bp - tp_bp
  base_precision <- safe_ratio(tp_bp, tp_bp + fp_bp)
  base_recall <- safe_ratio(tp_bp, tp_bp + fn_bp)
  base_rows[[label]] <- data.frame(
    length_group = label, length_bp = nominal_length,
    predicted_bp = predicted_bp, truth_bp = truth_bp,
    TP_bp = tp_bp, FP_bp = fp_bp, FN_bp = fn_bp,
    precision = base_precision, recall = base_recall,
    F1 = f1_value(base_precision, base_recall)
  )

  if (nrow(prediction) > 0 && any(!matched_prediction)) {
    tmp <- prediction[!matched_prediction, c("prediction_id", "chrom", "start", "end", "length_bp", "direction", "copy_number_A", "copy_number_B"), drop = FALSE]
    tmp$length_group <- label
    all_unmatched_predictions[[label]] <- tmp
  }
  if (any(!matched_truth)) {
    tmp <- truth[!matched_truth, c("truth_id", "chrom", "start", "end", "length_bp", "direction", "Sample_A_CN", "Sample_B_CN"), drop = FALSE]
    tmp$length_group <- label
    all_unmatched_truth[[label]] <- tmp
  }
}

event_metrics <- do.call(rbind, event_rows)
base_metrics <- do.call(rbind, base_rows)
rownames(event_metrics) <- NULL
rownames(base_metrics) <- NULL
write.csv(event_metrics, file.path(output_dir, "metrics_event_level.csv"), row.names = FALSE)
write.csv(base_metrics, file.path(output_dir, "metrics_base_level.csv"), row.names = FALSE)

combined <- rbind(
  data.frame(
    length_group = event_metrics$length_group,
    length_bp = event_metrics$length_bp,
    level = "Event level",
    precision = event_metrics$precision,
    recall = event_metrics$recall,
    F1 = event_metrics$F1
  ),
  data.frame(
    length_group = base_metrics$length_group,
    length_bp = base_metrics$length_bp,
    level = "Base-pair level",
    precision = base_metrics$precision,
    recall = base_metrics$recall,
    F1 = base_metrics$F1
  )
)
write.csv(combined, file.path(output_dir, "metrics_combined.csv"), row.names = FALSE)

if (length(all_matches) > 0) write.csv(do.call(rbind, all_matches), file.path(output_dir, "event_matches.csv"), row.names = FALSE)
if (length(all_unmatched_predictions) > 0) write.csv(do.call(rbind, all_unmatched_predictions), file.path(output_dir, "unmatched_predictions.csv"), row.names = FALSE)
if (length(all_unmatched_truth) > 0) write.csv(do.call(rbind, all_unmatched_truth), file.path(output_dir, "unmatched_truth.csv"), row.names = FALSE)

combined$length_group <- factor(combined$length_group, levels = labels)
combined$level <- factor(combined$level, levels = c("Event level", "Base-pair level"))
na_f1_rows <- combined[is.na(combined$F1), , drop = FALSE]
na_f1_rows$plot_y <- rep(0.04, nrow(na_f1_rows))
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}
f1_plot <- ggplot(combined, aes(x = length_group, y = F1, color = level, group = level)) +
  geom_hline(yintercept = seq(0, 1, 0.2), color = "#E5E7EB", linewidth = 0.30) +
  geom_line(linewidth = 0.75, na.rm = TRUE) +
  geom_point(size = 2.1, na.rm = TRUE) +
  geom_text(
    data = na_f1_rows,
    aes(y = plot_y, label = "NA"),
    family = "Arial", fontface = "bold", size = 2.1, show.legend = FALSE
  ) +
  scale_color_manual(values = c("Event level" = "#2C7FB8", "Base-pair level" = "#D95F0E")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "F1 score", color = NULL) +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    text = element_text(family = "Arial", size = 8),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 8),
    legend.key.width = unit(12, "pt"),
    legend.key.height = unit(8, "pt"),
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 10, margin = margin(r = 15)),
    axis.title.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = unit(3, "pt"),
    plot.margin = margin(t = 6, r = 6, b = 6, l = 10)
  )

grDevices::png(file.path(figure_dir, "F1_by_length.png"),
               width = 1650, height = 1050, res = 300, family = "Arial")
print(f1_plot)
grDevices::dev.off()
grDevices::cairo_pdf(file.path(figure_dir, "F1_by_length.pdf"),
                     width = 5.5, height = 3.5, family = "Arial")
print(f1_plot)
grDevices::dev.off()

method_lines <- c(
  "Coordinate convention: 0-based half-open intervals; interval length = end - start.",
  "Direction: dup when Sample_B CN - Sample_A CN > 0; del when the difference is < 0.",
  "Event candidate: same chromosome and direction, with prediction overlap fraction >= 0.80 OR truth overlap fraction >= 0.80.",
  "Event matching: greedy one-to-one, ranked by min(overlap fractions), then IoU, then overlap bp.",
  "Event metrics: TP = matched pairs; FP = unmatched predictions; FN = unmatched truth events.",
  "Base-pair metrics use only accepted one-to-one matches: TP_bp = sum of matched intersection bp; FP_bp = total predicted bp - TP_bp; FN_bp = total truth bp - TP_bp.",
  "Precision = TP/(TP+FP); recall = TP/(TP+FN); F1 = harmonic mean. The same formulas are applied to bp counts."
)
writeLines(method_lines, file.path(output_dir, "evaluation_method.txt"))

print(event_metrics)
print(base_metrics)

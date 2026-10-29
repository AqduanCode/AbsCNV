#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: call_cony_single_sample.R CONY_DIR DEPTH_SUM_TSV")

cony_dir <- normalizePath(args[[1]], mustWork = TRUE)
depth_path <- normalizePath(args[[2]], mustWork = TRUE)
work_dir <- file.path(cony_dir, "results")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

runtime_source <- file.path(work_dir, "CONY.R")
if (!file.copy(file.path(cony_dir, "software", "cony_compatible_runtime.R"), runtime_source,
               overwrite = TRUE)) stop("Unable to stage CONY.R")
source(runtime_source)

meta <- read.delim(file.path(cony_dir, "read_depth_10kb", "windows.10kb.meta.tsv"),
                   check.names = FALSE)
depth <- read.table(depth_path, header = FALSE,
                    col.names = c("chrom", "start0", "end0", "sum_A"))
if (nrow(meta) != nrow(depth) || !all(meta$chrom == depth$chrom) ||
    !all(meta$start0 == depth$start0) || !all(meta$end0 == depth$end0)) {
  stop("Window metadata and samtools bedcov output do not match")
}
dat <- cbind(meta, sum_A = depth$sum_A)

write_raw <- function(d) {
  denom <- d$nonAmb * d$width
  ard <- ifelse(denom > 0, d$sum_A / denom, NaN)
  raw <- data.frame(seq = d$chrom, start = d$start0 + 1L, end = d$end0,
                    width = d$width, nonAmb = d$nonAmb, GC = d$GC, ARD = ard)
  out <- sprintf("CONY.1-TempRegion.%s.Sample_A.SumUp.RD.txt", unique(d$chrom))
  write.table(raw, out, quote = FALSE, sep = " ", row.names = FALSE)
}

# CONY's upstream EstPar fixes the median read depth to CN2.  Keep that
# diploid baseline here; genome-wide mean absolute CN is not a baseline ploidy.
EstParPloidy <- function(TargetChr, SampleName = "Sample_A", NCN = 6L,
                        BaselinePloidy = 2) {
  path <- sprintf("CONY.3-TempRegion.%s.%s.SumUp.Single.UsedRD.txt",
                  TargetChr, SampleName)
  z <- read.table(path, header = TRUE)
  med0 <- median(na.omit(log(z$RD)))
  logv <- ((sort(log(z$RD), decreasing = TRUE)[round(nrow(z) * 0.025)] - med0) / 1.96)^2
  mean_p <- exp(med0 + logv / 2)
  var_p <- (exp(logv) - 1) * exp(2 * med0 + logv)
  scale <- (1:NCN) / BaselinePloidy
  om <- mean_p * scale
  ov <- var_p * scale
  gm <- log(om^2 / sqrt(ov + om^2))
  gv <- log(1 + ov / om^2)
  probs <- vapply(seq_len(NCN), function(i) dlnorm(z$RD, gm[i], sqrt(gv[i])),
                  numeric(nrow(z)))
  assigned <- max.col(probs, ties.method = "first")
  gp <- tabulate(assigned, nbins = NCN) / length(assigned)
  out <- data.frame(GroupCNV = 1:NCN, GroupMean = gm,
                    GroupVar = gv, GroupPro = gp)
  write.table(out, sprintf("CONY.3-GroupSumm.%s.%s.SumUp.Single.txt",
                           TargetChr, SampleName),
              quote = FALSE, sep = " ", row.names = FALSE)
}

chromosomes <- paste0("chr", 1:22)
all_rows <- list()
for (idx in seq_along(chromosomes)) {
  chr <- chromosomes[[idx]]
  final_file <- sprintf("CONY.Result.%s.Sample_A.SumUp.Single.CNV.txt", chr)
  if (!file.exists(final_file) || file.info(final_file)$size == 0) {
    message(sprintf("[%s] CONY baseline-ploidy-2 analysis started", chr))
    d <- dat[dat$chrom == chr, , drop = FALSE]
    if (nrow(d) == 0L) stop(sprintf("No windows found for %s", chr))
    write_raw(d)
    AdjRD(CRDMethod = "SumUp", TargetChr = chr, SampleName = "Sample_A")
    UsedRD(CRDMethod = "SumUp", AnaMethod = "Single", TargetChr = chr,
           SampleName = "Sample_A")
    EstParPloidy(TargetChr = chr, NCN = 6L, BaselinePloidy = 2)
    # Match the original baseline-3 run so baseline ploidy is the only
    # intentional model change between the two runs.
    set.seed(20261004L + idx)
    RunCONY(CRDMethod = "SumUp", AnaMethod = "Single", TargetChr = chr,
            SampleName = "Sample_A", RunTime = 300000, BurnN = 5000,
            RTN = 1000, BCPoint = 20, FragLength = 500000)
    ComResult(CRDMethod = "SumUp", AnaMethod = "Single", TargetChr = chr,
              SampleName = "Sample_A")
    message(sprintf("[%s] CONY baseline-ploidy-2 analysis complete", chr))
  }
  if (file.exists(final_file) && file.info(final_file)$size > 0) {
    z <- tryCatch(read.table(final_file, header = TRUE), error = function(e) NULL)
    if (!is.null(z) && nrow(z) > 0L) {
      z$chrom <- chr
      all_rows[[length(all_rows) + 1L]] <- z[, c("chrom", "start", "end", "CN")]
    }
  }
}

combined <- if (length(all_rows)) do.call(rbind, all_rows) else
  data.frame(chrom = character(), start = integer(), end = integer(), CN = integer())
write.table(combined, file.path(work_dir, "Sample_A.CONY.Single.10kb.BaselinePloidy2.CNV.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

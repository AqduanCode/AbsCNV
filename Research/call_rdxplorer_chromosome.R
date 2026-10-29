#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: call_rdxplorer_chromosome.R APP_HOME WORK_DIR CHROM")
}
apphome <- normalizePath(args[[1]], mustWork = TRUE)
wrkgdir <- normalizePath(args[[2]], mustWork = TRUE)
chrom <- args[[3]]
source(file.path(apphome, "rlib", "ewtMain.R"))
runEwtMain(
  apphome = apphome,
  wrkgdir = wrkgdir,
  countwinfile = paste0(chrom, "q20win100.txt"),
  chr = chrom,
  hg = "hg38",
  win = 100,
  gender = "F",
  baseCopy = 2,
  filter = 10,
  dirsep = "/"
)


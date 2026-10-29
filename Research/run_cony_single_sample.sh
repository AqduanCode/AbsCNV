#!/usr/bin/env bash
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
set -euo pipefail
if [[ $# -ne 3 ]]; then echo 'Usage: run_cony_single_sample.sh BAM REFERENCE_FASTA NEW_OUTPUT_DIR' >&2; exit 2; fi
bam=$1; ref=$2; out=$3
research=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
[[ -s "$bam" && -s "$bam.bai" && -s "$ref" ]]
mkdir -p "$out"/{software,read_depth_10kb,results,logs}
cp "$research/cony_compatible_runtime.R" "$out/software/cony_compatible_runtime.R"
python3 "$research/build_cony_depth_windows.py" "$ref" 10000 "$out/read_depth_10kb/windows.10kb.bed" "$out/read_depth_10kb/windows.10kb.meta.tsv"
samtools bedcov -Q 30 -j "$out/read_depth_10kb/windows.10kb.bed" "$bam" > "$out/read_depth_10kb/Sample_A.bedcov_Q30.tsv"
export CONY_WORKERS=${CONY_WORKERS:-48}
Rscript "$research/call_cony_single_sample.R" "$out" "$out/read_depth_10kb/Sample_A.bedcov_Q30.tsv" > "$out/logs/CONY.log" 2>&1
[[ -s "$out/results/Sample_A.CONY.Single.10kb.BaselinePloidy2.CNV.tsv" ]]

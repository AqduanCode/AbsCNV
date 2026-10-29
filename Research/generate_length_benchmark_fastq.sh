#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 LENGTH_LABEL" >&2
  exit 2
fi

LABEL=$1
PROJECT=${ABSCNV_PROJECT:-/home/data2/DAQ/HCC1395_like_CN_length_separate_benchmark_20260928}
EXP="$PROJECT/experiments/$LABEL"
INPUTS="$EXP/inputs"
FASTQ="$EXP/fastq"
CHROM_FASTQ="$EXP/cnvsim_output/chrom_fastq"
LOG_DIR="$EXP/logs"
QC="$EXP/qc"
REF_DIR=/home/data3/home/jjy/CNVsim/refs
SHARED_A=${ABSCNV_SAMPLE_A_FASTQ:-/home/data2/DAQ/HCC1395_like_CN_50events_benchmark_20260921/fastq}

mkdir -p "$FASTQ" "$CHROM_FASTQ" "$LOG_DIR" "$QC"
ln -sfn "$SHARED_A/Sample_A_R1.fastq.gz" "$FASTQ/Sample_A_R1.fastq.gz"
ln -sfn "$SHARED_A/Sample_A_R2.fastq.gz" "$FASTQ/Sample_A_R2.fastq.gz"

if [[ -s "$QC/simulation_completed_at.txt" && -s "$FASTQ/Sample_B_R1.fastq.gz" && -s "$FASTQ/Sample_B_R2.fastq.gz" ]]; then
  gzip -t "$FASTQ/Sample_B_R1.fastq.gz" "$FASTQ/Sample_B_R2.fastq.gz"
  echo "$LABEL simulation already complete"
  exit 0
fi

run_one() {
  local chrom=$1
  python3 "$RESEARCH_DIR/simulate_length_benchmark_reads.py" \
    Sample_B "$chrom" "$INPUTS/Sample_B.absolute_CN_segments.tsv" \
    "$REF_DIR/${chrom}.fa" "$CHROM_FASTQ" \
    "$LOG_DIR/Sample_B.${chrom}.segments.tsv" \
    > "$LOG_DIR/Sample_B.${chrom}.stdout.txt" 2>&1
}
export -f run_one
export PROJECT INPUTS REF_DIR CHROM_FASTQ LOG_DIR

printf 'chrom\n' > "$QC/simulation_jobs.tsv"
for n in $(seq 1 22); do printf 'chr%s\n' "$n" >> "$QC/simulation_jobs.tsv"; done
tail -n +2 "$QC/simulation_jobs.tsv" | xargs -n 1 -P 2 bash -c 'run_one "$0"'

: > "$FASTQ/Sample_B_R1.fastq.gz"
: > "$FASTQ/Sample_B_R2.fastq.gz"
for n in $(seq 1 22); do
  cat "$CHROM_FASTQ/Sample_B.chr${n}.R1.fastq.gz" >> "$FASTQ/Sample_B_R1.fastq.gz"
  cat "$CHROM_FASTQ/Sample_B.chr${n}.R2.fastq.gz" >> "$FASTQ/Sample_B_R2.fastq.gz"
done

gzip -t "$FASTQ/Sample_B_R1.fastq.gz" "$FASTQ/Sample_B_R2.fastq.gz"
sha256sum "$FASTQ/Sample_B_R1.fastq.gz" "$FASTQ/Sample_B_R2.fastq.gz" > "$QC/Sample_B.fastq.sha256"

python3 - "$LOG_DIR" "$QC/read_pair_summary.tsv" <<'PY'
from pathlib import Path
import csv, sys
log_dir, output = map(Path, sys.argv[1:])
total = 0
for path in log_dir.glob("Sample_B.chr*.segments.tsv"):
    with path.open() as handle:
        total += sum(int(row["read_pairs"]) for row in csv.DictReader(handle, delimiter="\t"))
output.write_text(f"sample\tread_pairs\nSample_B\t{total}\n")
PY

date '+%F %T %z' > "$QC/simulation_completed_at.txt"
find "$CHROM_FASTQ" -maxdepth 1 -type f -name 'Sample_B.chr*.fastq.gz' -delete
rmdir "$CHROM_FASTQ" 2>/dev/null || true
echo "$LABEL simulation complete"

#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

PROJECT=${ABSCNV_PROJECT:-/home/data2/DAQ/HCC1395_like_CN_50events_benchmark_20260921}
REF_DIR=/home/data3/home/jjy/CNVsim/refs
INPUTS="$PROJECT/inputs"
CHROM_FASTQ="$PROJECT/cnvsim_output/chrom_fastq"
LOG_DIR="$PROJECT/logs"
QC="$PROJECT/qc"
mkdir -p "$CHROM_FASTQ" "$LOG_DIR" "$QC" "$PROJECT/fastq"

python3 "$RESEARCH_DIR/build_paired_benchmark_truth.py" \
  "$PROJECT/source/cnv_gain_cn_median.txt" \
  "$PROJECT/source/cnv_loss_cn_median.txt" \
  "$REF_DIR" "$INPUTS"

cp "$INPUTS/Sample_B_vs_Sample_A.benchmark.tsv" "$PROJECT/benchmark/"
cp "$INPUTS/design_summary.tsv" "$PROJECT/benchmark/"

run_one() {
  local sample=$1 chrom=$2
  python3 "$RESEARCH_DIR/simulate_paired_benchmark_reads.py" \
    "$sample" "$chrom" "$INPUTS/${sample}.absolute_CN_segments.tsv" \
    "$REF_DIR/${chrom}.fa" "$CHROM_FASTQ" \
    "$LOG_DIR/${sample}.${chrom}.segments.tsv" \
    > "$LOG_DIR/${sample}.${chrom}.stdout.txt" 2>&1
}
export -f run_one
export PROJECT REF_DIR INPUTS CHROM_FASTQ LOG_DIR

printf 'sample\tchrom\n' > "$QC/jobs.tsv"
for sample in Sample_A Sample_B; do
  for n in $(seq 1 22); do
    printf '%s\tchr%s\n' "$sample" "$n" >> "$QC/jobs.tsv"
  done
done

tail -n +2 "$QC/jobs.tsv" | xargs -n 2 -P 8 bash -c 'run_one "$0" "$1"'

for sample in Sample_A Sample_B; do
  : > "$PROJECT/fastq/${sample}_R1.fastq.gz"
  : > "$PROJECT/fastq/${sample}_R2.fastq.gz"
  for n in $(seq 1 22); do
    cat "$CHROM_FASTQ/${sample}.chr${n}.R1.fastq.gz" >> "$PROJECT/fastq/${sample}_R1.fastq.gz"
    cat "$CHROM_FASTQ/${sample}.chr${n}.R2.fastq.gz" >> "$PROJECT/fastq/${sample}_R2.fastq.gz"
  done
done

gzip -t "$PROJECT"/fastq/*.fastq.gz
sha256sum "$PROJECT"/fastq/*.fastq.gz > "$QC/fastq.sha256"

python3 - "$LOG_DIR" "$QC/read_pair_summary.tsv" <<'PY'
from pathlib import Path
import csv, sys
log_dir, output = map(Path, sys.argv[1:])
with output.open('w') as out:
    out.write('sample\tread_pairs\n')
    for sample in ('Sample_A', 'Sample_B'):
        total = 0
        for path in log_dir.glob(f'{sample}.chr*.segments.tsv'):
            with path.open() as handle:
                total += sum(int(r['read_pairs']) for r in csv.DictReader(handle, delimiter='\t'))
        out.write(f'{sample}\t{total}\n')
PY

date '+%F %T %z' > "$QC/completed_at.txt"

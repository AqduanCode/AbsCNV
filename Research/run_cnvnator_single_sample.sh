#!/usr/bin/env bash
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
set -eo pipefail

source /home/anqiduan/miniconda3/etc/profile.d/conda.sh
conda activate cnvnator
set -u

PROJECT=${ABSCNV_PROJECT:-/home/data2/DAQ/HCC1395_like_Sample_A_single_sample_CNV_benchmark_20261004}
OUT="$PROJECT/Tools/CNVnator_1kb_fresh_baseline2"
BAM="$PROJECT/alignment/bam/Sample_A.sorted.bam"
REF_CHROMS="$PROJECT/reference/chromosomes"
ROOT_FILE="$OUT/Sample_A.fresh.root"
CALLS="$OUT/Sample_A.cnvnator.calls.tsv"
LOG_DIR="$OUT/logs"
STATUS="$PROJECT/PIPELINE_STATUS.tsv"
CNVNATOR="$CONDA_PREFIX/bin/cnvnator"
SAMTOOLS=/usr/local/bin/samtools
BIN=1000
CHROMS=()
for i in $(seq 1 22); do CHROMS+=("chr${i}"); done

mkdir -p "$OUT" "$LOG_DIR"

record() {
  printf '%s\t%s\t%s\t%s\n' "$(date '+%F %T %z')" "$1" "$2" "$3" >> "$STATUS"
}

[[ -x "$CNVNATOR" ]] || { record cnvnator_fresh_cn2 failed "missing executable"; exit 1; }
[[ -s "$BAM" && -s "$BAM.bai" ]] || { record cnvnator_fresh_cn2 failed "missing BAM or BAI"; exit 1; }
"$SAMTOOLS" quickcheck -v "$BAM"
for chrom in "${CHROMS[@]}"; do
  [[ -s "$REF_CHROMS/${chrom}.fa" ]] || {
    record cnvnator_fresh_cn2 failed "missing reference $chrom"
    exit 1
  }
done

record cnvnator_fresh_cn2 started "fresh tree/his/stat/partition/call; bin=1kb; absolute-CN baseline=2 in export"
rm -f "$ROOT_FILE.tmp" "$CALLS.tmp"

"$CNVNATOR" -root "$ROOT_FILE.tmp" -chrom "${CHROMS[@]}" -tree "$BAM" \
  > "$LOG_DIR/01_tree.log" 2>&1
[[ -s "$ROOT_FILE.tmp" ]] || { record cnvnator_fresh_cn2 failed "empty ROOT after tree"; exit 1; }

"$CNVNATOR" -root "$ROOT_FILE.tmp" -chrom "${CHROMS[@]}" -his "$BIN" -d "$REF_CHROMS" \
  > "$LOG_DIR/02_his.log" 2>&1
"$CNVNATOR" -root "$ROOT_FILE.tmp" -chrom "${CHROMS[@]}" -stat "$BIN" \
  > "$LOG_DIR/03_stat.log" 2>&1
"$CNVNATOR" -root "$ROOT_FILE.tmp" -chrom "${CHROMS[@]}" -partition "$BIN" \
  > "$LOG_DIR/04_partition.log" 2>&1
"$CNVNATOR" -root "$ROOT_FILE.tmp" -chrom "${CHROMS[@]}" -call "$BIN" \
  > "$CALLS.tmp" 2> "$LOG_DIR/05_call.stderr.log"

[[ -s "$CALLS.tmp" ]] || { record cnvnator_fresh_cn2 failed "empty call output"; exit 1; }
mv "$ROOT_FILE.tmp" "$ROOT_FILE"
mv "$CALLS.tmp" "$CALLS"
date '+%F %T %z' > "$OUT/COMPLETED.txt"
record cnvnator_fresh_cn2 complete "calls=$CALLS"

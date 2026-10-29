#!/usr/bin/env bash
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
set -eo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 LENGTH_LABEL" >&2
  exit 2
fi

source /home/anqiduan/miniconda3/etc/profile.d/conda.sh
conda activate cnvnator
set -u

LABEL=$1
PROJECT=${ABSCNV_PROJECT:-/home/data2/DAQ/HCC1395_like_CN_length_separate_benchmark_20260928}
EXP="$PROJECT/experiments/$LABEL"
WORK="$EXP/alignment_cnvnator_rd_1kb"
FASTQ="$EXP/fastq"
SHARED=${ABSCNV_SAMPLE_A_ALIGNMENT:-/home/data2/DAQ/HCC1395_like_CN_50events_benchmark_20260921/alignment_cnvnator_rd_1kb_20260921}
REF="$SHARED/reference/GRCh_chr1_22.fa"
REF_CHROMS="$SHARED/reference/chromosomes"
BAM_DIR="$WORK/bam"
ROOT_DIR="$WORK/cnvnator_root"
CSV_DIR="$WORK/rd_csv"
LOG_DIR="$WORK/logs"
TMP_DIR="$WORK/tmp"
QC_DIR="$WORK/qc"

BWA=/home/anqiduan/miniconda3/envs/hcc1395_wgs/bin/bwa-mem2
SAMTOOLS=/usr/local/bin/samtools
CNVNATOR="$CONDA_PREFIX/bin/cnvnator"
ROOT="$CONDA_PREFIX/bin/root"
BWA_THREADS=24
SORT_THREADS=12
SORT_MEM=2G
BIN=1000
CHROMS=()
for i in $(seq 1 22); do CHROMS+=("chr${i}"); done

mkdir -p "$BAM_DIR" "$ROOT_DIR" "$CSV_DIR" "$LOG_DIR" "$TMP_DIR" "$QC_DIR"
ln -sfn "$SHARED/reference" "$WORK/reference"
ln -sfn "$SHARED/cnvnator_root/Sample_A.root" "$ROOT_DIR/Sample_A.root"
ln -sfn "$SHARED/rd_csv/Sample_A.rd_1kb.csv" "$CSV_DIR/Sample_A.rd_1kb.csv"

if [[ -s "$QC_DIR/completed_at.txt" && -s "$ROOT_DIR/Sample_B.root" && -s "$CSV_DIR/Sample_B.rd_1kb.csv" ]]; then
  echo "$LABEL alignment/CNVnator already complete"
  exit 0
fi

available_bytes=$(df --output=avail -B1 "$EXP" | tail -n 1 | tr -d ' ')
if (( available_bytes < 120000000000 )); then
  echo "ERROR: less than 120 GB free before alignment: $available_bytes" >&2
  exit 1
fi

for exe in "$BWA" "$SAMTOOLS" "$CNVNATOR" "$ROOT"; do
  [[ -x "$exe" ]] || { echo "Missing executable: $exe" >&2; exit 1; }
done
[[ -s "$REF" && -s "$REF.bwt.2bit.64" ]] || { echo "Shared indexed reference missing" >&2; exit 1; }

R1="$FASTQ/Sample_B_R1.fastq.gz"
R2="$FASTQ/Sample_B_R2.fastq.gz"
BAM="$BAM_DIR/Sample_B.sorted.bam"
ROOT_FILE="$ROOT_DIR/Sample_B.root"
CSV="$CSV_DIR/Sample_B.rd_1kb.csv"
SAMPLE_TMP="$TMP_DIR/Sample_B"
RG='@RG\tID:Sample_B\tSM:Sample_B\tPL:ILLUMINA'
mkdir -p "$SAMPLE_TMP"

if [[ ! -s "$BAM" || ! -s "$BAM.bai" ]]; then
  "$BWA" mem -t "$BWA_THREADS" -R "$RG" "$REF" "$R1" "$R2" \
    2> "$LOG_DIR/Sample_B.bwa.log" \
    | "$SAMTOOLS" sort -@ "$SORT_THREADS" -m "$SORT_MEM" \
        -T "$SAMPLE_TMP/sort" -o "$BAM.tmp" - \
        2> "$LOG_DIR/Sample_B.sort.log"
  "$SAMTOOLS" quickcheck -v "$BAM.tmp"
  mv "$BAM.tmp" "$BAM"
  "$SAMTOOLS" index -@ "$SORT_THREADS" "$BAM"
fi

"$SAMTOOLS" quickcheck -v "$BAM"
"$SAMTOOLS" flagstat -@ "$SORT_THREADS" "$BAM" > "$QC_DIR/Sample_B.flagstat.txt"
"$SAMTOOLS" idxstats "$BAM" > "$QC_DIR/Sample_B.idxstats.tsv"

if [[ ! -e "$ROOT_DIR/Sample_B.tree.done" ]]; then
  rm -f "$ROOT_FILE"
  "$CNVNATOR" -root "$ROOT_FILE" -chrom "${CHROMS[@]}" -tree "$BAM" \
    > "$LOG_DIR/Sample_B.cnvnator_tree.log" 2>&1
  [[ -s "$ROOT_FILE" ]]
  touch "$ROOT_DIR/Sample_B.tree.done"
fi

if [[ ! -e "$ROOT_DIR/Sample_B.his_1000.done" ]]; then
  "$CNVNATOR" -root "$ROOT_FILE" -chrom "${CHROMS[@]}" -his "$BIN" -d "$REF_CHROMS" \
    > "$LOG_DIR/Sample_B.cnvnator_his_1000.log" 2>&1
  touch "$ROOT_DIR/Sample_B.his_1000.done"
fi

if [[ ! -s "$CSV" ]]; then
  rm -f "$CSV.tmp"
  "$ROOT" -l -b -q \
    "$RESEARCH_DIR/export_autosomal_read_depth.C(\"$ROOT_FILE\",\"$CSV.tmp\")" \
    > "$LOG_DIR/Sample_B.root_export.log" 2>&1
  [[ $(wc -l < "$CSV.tmp") -ge 2000000 ]]
  mv "$CSV.tmp" "$CSV"
fi

{
  printf 'sample\tbam_bytes\troot_bytes\tcsv_bytes\tcsv_rows\tbam_retained\n'
  printf 'Sample_B\t%s\t%s\t%s\t%s\tfalse\n' \
    "$(stat -c %s "$BAM")" "$(stat -c %s "$ROOT_FILE")" \
    "$(stat -c %s "$CSV")" "$(wc -l < "$CSV")"
} > "$QC_DIR/Sample_B.output_summary.tsv"

rm -f "$BAM" "$BAM.bai"
rm -rf "$SAMPLE_TMP"
date '+%F %T %z' > "$QC_DIR/completed_at.txt"
echo "$LABEL alignment/CNVnator/1kb RD complete"

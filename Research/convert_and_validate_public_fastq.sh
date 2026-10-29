#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

BASE=${ABSCNV_BASE:-/home/data2/DAQ/GSE118715_WGS_alignment_20260805}
BIN="$BASE/tools/sratoolkit.3.4.1-ubuntu64/bin"
SRA_DIR="$BASE/sra"
FASTQ_DIR="$BASE/fastq"
TMP_ROOT="$BASE/tmp/fasterq"
LOG_DIR="$BASE/logs/fasterq"
QC_DIR="$BASE/qc/fastq"
SUMMARY="$QC_DIR/fastq_validation_summary.tsv"
THREADS=32
MIN_FREE_GB=450

mkdir -p "$FASTQ_DIR" "$TMP_ROOT" "$LOG_DIR" "$QC_DIR"
printf 'run\texpected_spots\tr1_reads\tr2_reads\tr1_bytes_gz\tr2_bytes_gz\tr1_md5\tr2_md5\tgzip_test\n' > "$SUMMARY"

free_gb() {
  df -Pk "$BASE" | awk 'NR==2 {printf "%d\n", $4/1024/1024}'
}

convert_one() {
  local run=$1 expected_spots=$2
  local log="$LOG_DIR/${run}.log"
  local input="$SRA_DIR/$run/$run.sra"
  local r1="$FASTQ_DIR/${run}_1.fastq"
  local r2="$FASTQ_DIR/${run}_2.fastq"
  local lines1 lines2 reads1 reads2 first1 first2 last1 last2
  local bytes1 bytes2 md51 md52 available

  [[ -f "$QC_DIR/../sra_validate/${run}.VALIDATED" ]] || { echo "ERROR: SRA not validated: $run"; return 2; }
  [[ -f "$input" ]] || { echo "ERROR: SRA missing: $input"; return 2; }
  available=$(free_gb)
  (( available >= MIN_FREE_GB )) || { echo "ERROR: only ${available} GB free before $run"; return 75; }

  rm -rf -- "$TMP_ROOT/$run"
  mkdir -p "$TMP_ROOT/$run"
  printf '[%s] START fasterq run=%s expected_spots=%s free_gb=%s\n' "$(date '+%F %T %z')" "$run" "$expected_spots" "$available" >> "$log"

  "$BIN/fasterq-dump" --split-files --threads "$THREADS" --progress \
    --outdir "$FASTQ_DIR" --temp "$TMP_ROOT/$run" "$input" >> "$log" 2>&1

  [[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: FASTQ output missing for $run" >> "$log"; return 3; }
  lines1=$(wc -l < "$r1")
  lines2=$(wc -l < "$r2")
  (( lines1 % 4 == 0 && lines2 % 4 == 0 )) || { echo "ERROR: FASTQ line count not divisible by 4" >> "$log"; return 4; }
  reads1=$((lines1 / 4))
  reads2=$((lines2 / 4))
  [[ "$reads1" == "$expected_spots" && "$reads2" == "$expected_spots" ]] || {
    printf 'ERROR: read count mismatch expected=%s r1=%s r2=%s\n' "$expected_spots" "$reads1" "$reads2" >> "$log"
    return 5
  }

  first1=$(head -1 "$r1" | awk '{print $1}')
  first2=$(head -1 "$r2" | awk '{print $1}')
  last1=$(tail -4 "$r1" | head -1 | awk '{print $1}')
  last2=$(tail -4 "$r2" | head -1 | awk '{print $1}')
  [[ "$first1" == "$first2" && "$last1" == "$last2" ]] || {
    printf 'ERROR: paired read IDs disagree first=(%s,%s) last=(%s,%s)\n' "$first1" "$first2" "$last1" "$last2" >> "$log"
    return 6
  }

  printf '[%s] read_count_and_pair_boundary_ids=PASS; compressing\n' "$(date '+%F %T %z')" >> "$log"
  pigz -p 24 "$r1" & c1=$!
  pigz -p 24 "$r2" & c2=$!
  compression_failed=0
  wait "$c1" || compression_failed=1
  wait "$c2" || compression_failed=1
  (( compression_failed == 0 )) || { echo 'ERROR: pigz failed' >> "$log"; return 7; }

  pigz -t "$r1.gz"
  pigz -t "$r2.gz"
  bytes1=$(stat -c %s "$r1.gz")
  bytes2=$(stat -c %s "$r2.gz")
  md51=$(md5sum "$r1.gz" | awk '{print $1}')
  md52=$(md5sum "$r2.gz" | awk '{print $1}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "$run" "$expected_spots" "$reads1" "$reads2" "$bytes1" "$bytes2" "$md51" "$md52" >> "$SUMMARY"
  touch "$QC_DIR/${run}.FASTQ_GZ_VALIDATED"
  rm -rf -- "$TMP_ROOT/$run"
  printf '[%s] COMPLETE run=%s fastq_gz_validation=PASS\n' "$(date '+%F %T %z')" "$run" >> "$log"
}

# Sequential conversion limits peak disk use.
convert_one SRR7707725 431609169
convert_one SRR7707727 459477815
convert_one SRR7707729 449005612
printf '[%s] RESULT=ALL_FASTQ_GZ_VALIDATED\n' "$(date '+%F %T %z')"

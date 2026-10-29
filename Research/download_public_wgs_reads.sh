#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

BASE=${ABSCNV_BASE:-/home/data2/DAQ/GSE118715_WGS_alignment_20260805}
BIN="$BASE/tools/sratoolkit.3.4.1-ubuntu64/bin"
SRA_DIR="$BASE/sra"
LOG_DIR="$BASE/logs/prefetch"
QC_DIR="$BASE/qc/sra_validate"

mkdir -p "$SRA_DIR" "$LOG_DIR" "$QC_DIR"

download_validate() {
  local run=$1 expected_size=$2 expected_md5=$3
  local log="$LOG_DIR/${run}.log"
  local sra_file

  printf '[%s] START run=%s\n' "$(date '+%F %T %z')" "$run" >> "$log"
  "$BIN/prefetch" --max-size u --progress --output-directory "$SRA_DIR" "$run" >> "$log" 2>&1

  sra_file="$SRA_DIR/$run/$run.sra"
  if [[ ! -f "$sra_file" ]]; then
    sra_file="$SRA_DIR/$run/$run"
  fi
  [[ -f "$sra_file" ]] || { printf '[%s] ERROR missing SRA file\n' "$(date '+%F %T %z')" >> "$log"; return 2; }

  local actual_size actual_md5
  actual_size=$(stat -c %s "$sra_file")
  actual_md5=$(md5sum "$sra_file" | awk '{print $1}')
  printf 'expected_size=%s actual_size=%s\n' "$expected_size" "$actual_size" | tee "$QC_DIR/${run}.checks.txt" >> "$log"
  printf 'expected_md5=%s actual_md5=%s\n' "$expected_md5" "$actual_md5" | tee -a "$QC_DIR/${run}.checks.txt" >> "$log"
  [[ "$actual_size" == "$expected_size" ]] || { printf '[%s] ERROR size mismatch\n' "$(date '+%F %T %z')" >> "$log"; return 3; }
  [[ "$actual_md5" == "$expected_md5" ]] || { printf '[%s] ERROR md5 mismatch\n' "$(date '+%F %T %z')" >> "$log"; return 4; }

  "$BIN/vdb-validate" "$SRA_DIR/$run" > "$QC_DIR/${run}.vdb-validate.txt" 2>&1
  touch "$QC_DIR/${run}.VALIDATED"
  printf '[%s] COMPLETE run=%s size_md5_vdb=PASS\n' "$(date '+%F %T %z')" "$run" >> "$log"
}

download_validate SRR7707725 64751574518 18d20e53d2cdbfa013246462c6ef397b & p1=$!
download_validate SRR7707727 68997701879 9f4c7684f10f631848fc1b5bfb048512 & p2=$!
download_validate SRR7707729 66287565971 1f6146a8ad1f724815daee5cbba9b67e & p3=$!

failed=0
wait "$p1" || failed=1
wait "$p2" || failed=1
wait "$p3" || failed=1
if (( failed != 0 )); then
  printf '[%s] RESULT=FAIL\n' "$(date '+%F %T %z')"
  exit 1
fi
printf '[%s] RESULT=ALL_VALIDATED\n' "$(date '+%F %T %z')"

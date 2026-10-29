#!/usr/bin/env bash
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
set -Eeo pipefail
source /home/anqiduan/miniconda3/etc/profile.d/conda.sh
conda activate cnvnator
set -u

BASE=${ABSCNV_BASE:-/home/data2/DAQ/GSE118715_WGS_alignment_20260805}
OUT="$BASE/CNV_results/CNVnator_1kb"
BAM_DIR="$BASE/bam"
: "${ABSCNV_REFERENCE_CHROMS:?Set ABSCNV_REFERENCE_CHROMS to the matching chr1.fa through chr22.fa directory}"
REF_SOURCE="$ABSCNV_REFERENCE_CHROMS"
CNVNATOR="$CONDA_PREFIX/bin/cnvnator"
ROOT="$CONDA_PREFIX/bin/root"
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
MACRO="$RESEARCH_DIR/export_autosomal_read_depth.C"
BIN=1000
CHROMS=()
for i in $(seq 1 22); do CHROMS+=("chr${i}"); done

mkdir -p "$OUT"/{root,rd_csv,logs,qc,tmp,reference/chromosomes}
for chrom in "${CHROMS[@]}"; do
  ln -sfn "$REF_SOURCE/${chrom}.fa" "$OUT/reference/chromosomes/${chrom}.fa"
done

run_one() {
  local sample=$1 run=$2
  local bam="$BAM_DIR/${sample}.${run}.hg38.sorted.bam"
  local root_file="$OUT/root/${sample}.root"
  local csv="$OUT/rd_csv/${sample}.rd_1kb.csv"
  local log="$OUT/logs/${sample}.log"
  exec > >(tee -a "$log") 2>&1
  echo "[$(date '+%F %T')] START $sample"
  /home/anqiduan/miniconda3/envs/hcc1395_wgs/bin/samtools quickcheck -v "$bam"
  if [[ ! -f "$OUT/qc/${sample}.tree.done" ]]; then
    rm -f "$root_file"
    "$CNVNATOR" -root "$root_file" -chrom "${CHROMS[@]}" -tree "$bam"
    [[ -s "$root_file" ]]
    touch "$OUT/qc/${sample}.tree.done"
  fi
  if [[ ! -f "$OUT/qc/${sample}.his_${BIN}.done" ]]; then
    "$CNVNATOR" -root "$root_file" -chrom "${CHROMS[@]}" -his "$BIN" -d "$OUT/reference/chromosomes"
    touch "$OUT/qc/${sample}.his_${BIN}.done"
  fi
  if [[ ! -s "$csv" ]]; then
    rm -f "$csv.tmp"
    "$ROOT" -l -b -q "$MACRO(\"$root_file\",\"$csv.tmp\")"
    rows=$(wc -l < "$csv.tmp")
    (( rows > 2000000 )) || { echo "ERROR: only $rows CSV rows"; return 1; }
    mv "$csv.tmp" "$csv"
  fi
  printf 'sample\tbam_bytes\troot_bytes\tcsv_bytes\tcsv_rows\n%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$(stat -c %s "$bam")" "$(stat -c %s "$root_file")" \
    "$(stat -c %s "$csv")" "$(wc -l < "$csv")" > "$OUT/qc/${sample}.summary.tsv"
  touch "$OUT/qc/${sample}.CNVNATOR_1KB_VALIDATED"
  echo "[$(date '+%F %T')] PASS $sample"
}

run_one MCF7 SRR7707725 & p1=$!
run_one TAMR SRR7707727 & p2=$!
run_one FASR SRR7707729 & p3=$!
rc=0
wait "$p1" || rc=1
wait "$p2" || rc=1
wait "$p3" || rc=1
echo "[$(date '+%F %T')] CNVNATOR_WORKERS_DONE rc=$rc"
(( rc == 0 )) && touch "$OUT/COMPLETED.txt"
exit "$rc"

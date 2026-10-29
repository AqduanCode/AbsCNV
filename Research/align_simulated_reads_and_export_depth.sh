#!/usr/bin/env bash
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
set -eo pipefail

# ROOT/Cling needs the compiler and C++ standard-library headers supplied by
# the cnvnator Conda environment, not only the root executable from that env.
source /home/anqiduan/miniconda3/etc/profile.d/conda.sh
conda activate cnvnator
set -u

BASE=${ABSCNV_BASE:-/home/data2/DAQ/HCC1395_like_CN_50events_benchmark_20260921}
WORK="$BASE/alignment_cnvnator_rd_1kb_20260921"
FASTQ="$BASE/fastq"
REF_SOURCE=/home/data3/home/jjy/CNVsim/refs
REF_DIR="$WORK/reference"
REF="$REF_DIR/GRCh_chr1_22.fa"
BAM_DIR="$WORK/bam"
ROOT_DIR="$WORK/cnvnator_root"
CSV_DIR="$WORK/rd_csv"
SCRIPT_DIR="$WORK/scripts"
LOG_DIR="$WORK/logs"
TMP_DIR="$WORK/tmp"
QC_DIR="$WORK/qc"

BWA=/home/anqiduan/miniconda3/envs/hcc1395_wgs/bin/bwa-mem2
SAMTOOLS=/usr/local/bin/samtools
CNVNATOR="$CONDA_PREFIX/bin/cnvnator"
ROOT="$CONDA_PREFIX/bin/root"

BWA_THREADS=48
SORT_THREADS=16
SORT_MEM=2G
BIN=1000
CHROMS=()
for i in $(seq 1 22); do CHROMS+=("chr${i}"); done

mkdir -p "$REF_DIR/chromosomes" "$BAM_DIR" "$ROOT_DIR" "$CSV_DIR" \
  "$SCRIPT_DIR" "$LOG_DIR" "$TMP_DIR" "$QC_DIR"

log() {
  printf '[%s] %s\n' "$(date '+%F %T %z')" "$*"
}

for exe in "$BWA" "$SAMTOOLS" "$CNVNATOR" "$ROOT"; do
  if [[ ! -x "$exe" ]]; then
    log "ERROR: executable not found: $exe"
    exit 1
  fi
done

for sample in Sample_A Sample_B; do
  for mate in 1 2; do
    fq="$FASTQ/${sample}_R${mate}.fastq.gz"
    if [[ ! -s "$fq" ]]; then
      log "ERROR: FASTQ missing or empty: $fq"
      exit 1
    fi
  done
done

if [[ ! -s "$REF" ]]; then
  log "Building chr1-chr22 combined reference: $REF"
  tmp_ref="$REF.tmp.$$"
  : > "$tmp_ref"
  for chrom in "${CHROMS[@]}"; do
    src="$REF_SOURCE/${chrom}.fa"
    if [[ ! -s "$src" ]]; then
      log "ERROR: missing reference chromosome: $src"
      rm -f "$tmp_ref"
      exit 1
    fi
    cat "$src" >> "$tmp_ref"
    ln -sfn "$src" "$REF_DIR/chromosomes/${chrom}.fa"
  done
  mv "$tmp_ref" "$REF"
fi

for chrom in "${CHROMS[@]}"; do
  ln -sfn "$REF_SOURCE/${chrom}.fa" "$REF_DIR/chromosomes/${chrom}.fa"
done

if [[ ! -s "$REF.0123" || ! -s "$REF.bwt.2bit.64" ]]; then
  log "Indexing reference with bwa-mem2"
  "$BWA" index "$REF" > "$LOG_DIR/bwa_index.log" 2>&1
fi

if [[ ! -s "$REF.fai" ]]; then
  log "Creating FASTA index"
  "$SAMTOOLS" faidx "$REF"
fi

"$SAMTOOLS" dict -o "$REF_DIR/GRCh_chr1_22.dict.tmp" "$REF" 2>/dev/null || true
if [[ -s "$REF_DIR/GRCh_chr1_22.dict.tmp" ]]; then
  mv "$REF_DIR/GRCh_chr1_22.dict.tmp" "$REF_DIR/GRCh_chr1_22.dict"
else
  rm -f "$REF_DIR/GRCh_chr1_22.dict.tmp"
fi

run_sample() {
  local sample="$1"
  local r1="$FASTQ/${sample}_R1.fastq.gz"
  local r2="$FASTQ/${sample}_R2.fastq.gz"
  local bam="$BAM_DIR/${sample}.sorted.bam"
  local bai="$bam.bai"
  local root_file="$ROOT_DIR/${sample}.root"
  local csv="$CSV_DIR/${sample}.rd_1kb.csv"
  local sample_tmp="$TMP_DIR/$sample"
  local rg="@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA"

  mkdir -p "$sample_tmp"

  if [[ ! -s "$bam" || ! -s "$bai" ]]; then
    log "$sample: aligning paired FASTQs and coordinate-sorting BAM"
    rm -f "$bam.tmp" "$bam.tmp.bai"
    "$BWA" mem -t "$BWA_THREADS" -R "$rg" "$REF" "$r1" "$r2" \
      2> "$LOG_DIR/${sample}.bwa.log" \
      | "$SAMTOOLS" sort -@ "$SORT_THREADS" -m "$SORT_MEM" \
          -T "$sample_tmp/sort" -o "$bam.tmp" - \
          2> "$LOG_DIR/${sample}.sort.log"
    "$SAMTOOLS" quickcheck -v "$bam.tmp"
    mv "$bam.tmp" "$bam"
    "$SAMTOOLS" index -@ "$SORT_THREADS" "$bam"
  else
    log "$sample: existing BAM and index found; skipping alignment"
  fi

  "$SAMTOOLS" quickcheck -v "$bam"
  if [[ ! -s "$QC_DIR/${sample}.flagstat.txt" ]]; then
    "$SAMTOOLS" flagstat -@ "$SORT_THREADS" "$bam" > "$QC_DIR/${sample}.flagstat.txt"
  fi
  "$SAMTOOLS" idxstats "$bam" > "$QC_DIR/${sample}.idxstats.tsv"

  if [[ ! -e "$ROOT_DIR/${sample}.tree.done" ]]; then
    log "$sample: importing BAM read-depth tree into CNVnator ROOT"
    rm -f "$root_file"
    "$CNVNATOR" -root "$root_file" -chrom "${CHROMS[@]}" -tree "$bam" \
      > "$LOG_DIR/${sample}.cnvnator_tree.log" 2>&1
    [[ -s "$root_file" ]]
    touch "$ROOT_DIR/${sample}.tree.done"
  fi

  if [[ ! -e "$ROOT_DIR/${sample}.his_${BIN}.done" ]]; then
    log "$sample: calculating CNVnator RD histograms at ${BIN}-bp bins"
    "$CNVNATOR" -root "$root_file" -chrom "${CHROMS[@]}" \
      -his "$BIN" -d "$REF_DIR/chromosomes" \
      > "$LOG_DIR/${sample}.cnvnator_his_${BIN}.log" 2>&1
    touch "$ROOT_DIR/${sample}.his_${BIN}.done"
  fi

  if [[ ! -s "$csv" ]]; then
    log "$sample: exporting ${BIN}-bp RD histograms to CSV"
    rm -f "$csv.tmp"
    "$ROOT" -l -b -q \
      "$SCRIPT_DIR/export_autosomal_read_depth.C(\"$root_file\",\"$csv.tmp\")" \
      > "$LOG_DIR/${sample}.root_export.log" 2>&1
    if [[ $(wc -l < "$csv.tmp") -lt 2000000 ]]; then
      log "ERROR: exported CSV has unexpectedly few rows: $csv.tmp"
      exit 1
    fi
    mv "$csv.tmp" "$csv"
  fi

  {
    printf 'sample\tbam_bytes\troot_bytes\tcsv_bytes\tcsv_rows\n'
    printf '%s\t%s\t%s\t%s\t%s\n' "$sample" \
      "$(stat -c %s "$bam")" "$(stat -c %s "$root_file")" \
      "$(stat -c %s "$csv")" "$(wc -l < "$csv")"
  } > "$QC_DIR/${sample}.output_summary.tsv"

  rm -rf "$sample_tmp"
  log "$sample: complete"
}

log "Pipeline started"
run_sample Sample_A
run_sample Sample_B
date '+%F %T %z' > "$QC_DIR/completed_at.txt"
log "Pipeline completed"

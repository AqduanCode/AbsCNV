#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

PROJECT=${ABSCNV_PROJECT:-/home/data2/DAQ/HCC1395_like_Sample_A_single_sample_CNV_benchmark_20261004}
SOURCE=/home/data2/DAQ/HCC1395_like_CN_50events_benchmark_20260921
REF_SOURCE=/home/data3/home/jjy/CNVsim/refs
FASTQ="$SOURCE/fastq"
REF_DIR="$PROJECT/reference"
REF="$REF_DIR/GRCh_chr1_22.fa"
BAM_DIR="$PROJECT/alignment/bam"
LOG_DIR="$PROJECT/alignment/logs"
TMP_DIR="$PROJECT/alignment/tmp/Sample_A"
QC_DIR="$PROJECT/alignment/qc"
STATUS="$PROJECT/PIPELINE_STATUS.tsv"

BWA=/home/anqiduan/miniconda3/envs/hcc1395_wgs/bin/bwa-mem2
SAMTOOLS=/usr/local/bin/samtools
BWA_THREADS=48
SORT_THREADS=16
SORT_MEM=2G

mkdir -p "$REF_DIR/chromosomes" "$BAM_DIR" "$LOG_DIR" "$TMP_DIR" "$QC_DIR"

record() {
  printf '%s\t%s\t%s\t%s\n' "$(date '+%F %T %z')" "$1" "$2" "$3" >> "$STATUS"
}

if [[ ! -s "$STATUS" ]]; then
  printf 'time\tstage\tstatus\tdetail\n' > "$STATUS"
fi

for exe in "$BWA" "$SAMTOOLS"; do
  [[ -x "$exe" ]] || { record alignment failed "missing executable: $exe"; exit 1; }
done
for mate in 1 2; do
  fq="$FASTQ/Sample_A_R${mate}.fastq.gz"
  [[ -s "$fq" ]] || { record alignment failed "missing FASTQ: $fq"; exit 1; }
done

if [[ ! -s "$REF" ]]; then
  record reference started "concatenate GRCh38 chr1-chr22 from $REF_SOURCE"
  tmp_ref="$REF.tmp.$$"
  : > "$tmp_ref"
  for i in $(seq 1 22); do
    chrom="chr${i}"
    src="$REF_SOURCE/${chrom}.fa"
    [[ -s "$src" ]] || { record reference failed "missing chromosome FASTA: $src"; exit 1; }
    cat "$src" >> "$tmp_ref"
    ln -sfn "$src" "$REF_DIR/chromosomes/${chrom}.fa"
  done
  mv "$tmp_ref" "$REF"
  record reference complete "combined reference created"
fi

for i in $(seq 1 22); do
  chrom="chr${i}"
  ln -sfn "$REF_SOURCE/${chrom}.fa" "$REF_DIR/chromosomes/${chrom}.fa"
done

if [[ ! -s "$REF.0123" || ! -s "$REF.bwt.2bit.64" ]]; then
  record reference_index started "bwa-mem2 index"
  "$BWA" index "$REF" > "$LOG_DIR/bwa_index.log" 2>&1
  record reference_index complete "bwa-mem2 index complete"
fi
if [[ ! -s "$REF.fai" ]]; then
  "$SAMTOOLS" faidx "$REF"
fi

bam="$BAM_DIR/Sample_A.sorted.bam"
r1="$FASTQ/Sample_A_R1.fastq.gz"
r2="$FASTQ/Sample_A_R2.fastq.gz"
rg='@RG\tID:Sample_A\tSM:Sample_A\tPL:ILLUMINA'

if [[ ! -s "$bam" || ! -s "$bam.bai" ]]; then
  record alignment started "bwa-mem2 48 threads; samtools coordinate sort 16 threads"
  rm -f "$bam.tmp" "$bam.tmp.bai"
  "$BWA" mem -t "$BWA_THREADS" -R "$rg" "$REF" "$r1" "$r2" \
    2> "$LOG_DIR/Sample_A.bwa.log" \
    | "$SAMTOOLS" sort -@ "$SORT_THREADS" -m "$SORT_MEM" \
        -T "$TMP_DIR/sort" -o "$bam.tmp" - \
        2> "$LOG_DIR/Sample_A.sort.log"
  "$SAMTOOLS" quickcheck -v "$bam.tmp"
  mv "$bam.tmp" "$bam"
  "$SAMTOOLS" index -@ "$SORT_THREADS" "$bam"
  record alignment complete "BAM and BAI created"
else
  record alignment skipped "existing BAM and BAI retained"
fi

"$SAMTOOLS" quickcheck -v "$bam"
"$SAMTOOLS" flagstat -@ "$SORT_THREADS" "$bam" > "$QC_DIR/Sample_A.flagstat.txt"
"$SAMTOOLS" idxstats "$bam" > "$QC_DIR/Sample_A.idxstats.tsv"
sha256sum "$bam.bai" > "$QC_DIR/Sample_A.bai.sha256"
stat -c 'bam_bytes\t%s' "$bam" > "$QC_DIR/Sample_A.output_summary.tsv"
printf 'completed_at\t%s\n' "$(date '+%F %T %z')" >> "$QC_DIR/Sample_A.output_summary.tsv"
touch "$QC_DIR/ALIGNMENT_COMPLETE"
rm -rf "$TMP_DIR"

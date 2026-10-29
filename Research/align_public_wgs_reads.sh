#!/usr/bin/env bash
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
set -Eeuo pipefail

BASE=${ABSCNV_BASE:-/home/data2/DAQ/GSE118715_WGS_alignment_20260805}
FASTQ="$BASE/fastq"
BAM="$BASE/bam"
QC="$BASE/qc/alignment"
LOG="$BASE/logs/alignment"
TMP="$BASE/tmp/alignment"
REF=/home/data2/DAQ/reference/human/hg38.fa
BWA=/home/anqiduan/miniconda3/envs/hcc1395_wgs/bin/bwa-mem2
SAMTOOLS=/home/anqiduan/miniconda3/envs/hcc1395_wgs/bin/samtools
THREADS=16
MIN_FREE_GIB=120

mkdir -p "$BAM" "$QC" "$LOG" "$TMP"

align_one() {
    local sample=$1 run=$2 expected_spots=$3
    local r1="$FASTQ/${run}_1.fastq.gz"
    local r2="$FASTQ/${run}_2.fastq.gz"
    local marker="$BASE/qc/fastq/${run}.FASTQ_GZ_VALIDATED"
    local out="$BAM/${sample}.${run}.hg38.sorted.bam"
    local part="${out}.part"
    local sample_log="$LOG/${sample}.${run}.align.log"
    local done_marker="$QC/${sample}.${run}.ALIGNMENT_VALIDATED"

    exec > >(tee -a "$sample_log") 2>&1
    echo "[$(date '+%F %T')] waiting for validated FASTQ: $run"
    while [[ ! -f "$marker" ]]; do
        sleep 60
    done
    echo "[$(date '+%F %T')] FASTQ marker found: $marker"

    if [[ -f "$done_marker" ]] && "$SAMTOOLS" quickcheck -v "$out"; then
        echo "[$(date '+%F %T')] already validated, skipping: $out"
        return 0
    fi

    [[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: paired FASTQ missing"; return 1; }
    gzip -t "$r1" "$r2"

    local free_gib
    free_gib=$(df -BG --output=avail "$BASE" | tail -1 | tr -dc '0-9')
    if (( free_gib < MIN_FREE_GIB )); then
        echo "ERROR: only ${free_gib} GiB free; threshold is ${MIN_FREE_GIB} GiB"
        return 1
    fi

    rm -f "$part" "${part}.bai"
    echo "[$(date '+%F %T')] START $sample $run; free=${free_gib}GiB"
    "$BWA" mem -t "$THREADS" \
        -R "@RG\\tID:${run}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:${run}" \
        "$REF" "$r1" "$r2" \
      | "$SAMTOOLS" sort -@ "$THREADS" -m 2G \
          -T "$TMP/${sample}.${run}" -o "$part" -

    "$SAMTOOLS" quickcheck -v "$part"
    mv "$part" "$out"
    "$SAMTOOLS" index -@ "$THREADS" "$out"
    "$SAMTOOLS" flagstat -@ "$THREADS" "$out" > "$QC/${sample}.${run}.flagstat.txt"
    "$SAMTOOLS" idxstats "$out" > "$QC/${sample}.${run}.idxstats.txt"

    local expected_primary actual_primary
    expected_primary=$((expected_spots * 2))
    actual_primary=$(awk '$4=="primary" {print $1; exit}' "$QC/${sample}.${run}.flagstat.txt")
    if [[ "$actual_primary" != "$expected_primary" ]]; then
        echo "ERROR: primary records=$actual_primary, expected=$expected_primary"
        return 1
    fi

    touch "$done_marker"
    echo "[$(date '+%F %T')] PASS $sample $run; primary=$actual_primary; bam=$(stat -c %s "$out")"
}

align_one MCF7 SRR7707725 431609169 & p1=$!
align_one TAMR SRR7707727 459477815 & p2=$!
align_one FASR SRR7707729 449005612 & p3=$!

rc=0
wait "$p1" || rc=1
wait "$p2" || rc=1
wait "$p3" || rc=1
echo "[$(date '+%F %T')] ALL_ALIGNMENT_WORKERS_FINISHED rc=$rc"
exit "$rc"

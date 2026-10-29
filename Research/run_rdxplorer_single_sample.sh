#!/usr/bin/env bash
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
set -euo pipefail
if [[ $# -ne 4 ]]; then echo 'Usage: run_rdxplorer_single_sample.sh BAM REF_FASTA PATCHED_APP NEW_OUTPUT_DIR' >&2; exit 2; fi
bam=$1; ref=$2; app=$3; out=$4
research=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
[[ -s "$bam" && -s "$bam.bai" && -s "$ref.fai" && -d "$app/rlib/gc100_hg38" && -d "$app/rlib/gap_hg38" ]]
mkdir -p "$out/results" "$out/logs"
for n in $(seq 1 22); do
  chrom="chr$n"
  length=$(awk -v c="$chrom" '$1==c {print $2}' "$ref.fai")
  [[ -n "$length" ]]
  bins=$((length / 100))
  samtools depth -aa -Q 20 -q 0 -r "$chrom" "$bam" \
    | awk -v N="$bins" 'BEGIN{b=0;s=0} {idx=int(($2-1)/100); while(b<idx && b<N){print s; b++; s=0} if(idx<N)s+=$3} END{while(b<N){print s; b++; s=0}}' \
    > "$out/results/${chrom}q20win100.txt"
  [[ $(wc -l < "$out/results/${chrom}q20win100.txt") -eq "$bins" ]]
  Rscript "$research/call_rdxplorer_chromosome.R" "$app" "$out/results" "$chrom" > "$out/logs/$chrom.log" 2>&1
  [[ -s "$out/results/$chrom.sum.filt" ]]
done

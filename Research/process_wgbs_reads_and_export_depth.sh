#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Reconstructed from the manuscript and recovered CNVnator workflow, 2026-10-28.
# This is a new configurable entry point, not a recovered historical master script.
# The manuscript reports fastp 0.23.1 / Bismark 0.22.3, whereas two inspected
# sorted U251 BAM headers record Bismark 0.24.2. See preprocessing provenance.
if [[ $# -ne 6 ]]; then
  echo 'Usage: process_wgbs_reads_and_export_depth.sh SAMPLE R1_FASTQ R2_FASTQ BISMARK_GENOME SPLIT_FASTA_DIR OUTPUT_DIR' >&2
  exit 2
fi
sample=$1; r1=$2; r2=$3; genome=$4; split_fasta=$5; out=$6
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$out/clean" "$out/alignment" "$out/deduplicated" "$out/methylation" "$out/rd" "$out/logs"
fastp --version > "$out/logs/fastp_version.txt" 2>&1
bismark --version > "$out/logs/bismark_version.txt" 2>&1
fastp -i "$r1" -I "$r2" -o "$out/clean/${sample}_1.clean.fq.gz" -O "$out/clean/${sample}_2.clean.fq.gz" \
  -q 5 -u 50 -f 10 -F 10 -n 15 -l 140 --overlap_diff_limit 1 --overlap_diff_percent_limit 10 \
  --json "$out/logs/fastp.json" --html "$out/logs/fastp.html" > "$out/logs/fastp.log" 2>&1
bismark --genome "$genome" -1 "$out/clean/${sample}_1.clean.fq.gz" -2 "$out/clean/${sample}_2.clean.fq.gz" \
  --output_dir "$out/alignment" > "$out/logs/bismark.log" 2>&1
shopt -s nullglob
aligned=("$out/alignment/"*_bismark_bt2_pe.bam)
[[ ${#aligned[@]} -eq 1 ]] || { echo 'Expected exactly one paired-end Bismark BAM' >&2; exit 1; }
deduplicate_bismark --paired --bam "${aligned[0]}" --output_dir "$out/deduplicated" > "$out/logs/deduplicate_bismark.log" 2>&1
dedup=("$out/deduplicated/"*.deduplicated.bam)
[[ ${#dedup[@]} -eq 1 ]] || { echo 'Expected exactly one deduplicated BAM' >&2; exit 1; }
root_file="$out/rd/${sample}.root"
[[ ! -e "$root_file" ]] || { echo 'Use a new output directory; ROOT already exists' >&2; exit 1; }
cnvnator -root "$root_file" -tree "${dedup[0]}"
cnvnator -root "$root_file" -his 1000 -d "$split_fasta"
cnvnator -root "$root_file" -stat 1000
cnvnator -root "$root_file" -partition 1000
cnvnator -root "$root_file" -call 1000 > "$out/rd/${sample}.1000-cnvout.txt"
root -l -b -q "$script_dir/export_gc_corrected_read_depth.C(\"$root_file\",\"$out/rd/${sample}.csv\")"
[[ -s "$out/rd/${sample}.csv" ]]
if [[ ${RUN_METHYLATION:-0} == 1 ]]; then
  bismark_methylation_extractor -p --gzip --comprehensive --no_overlap --bedGraph --counts \
    --report --cytosine_report --multicore 3 --genome_folder "$genome" \
    --output "$out/methylation" "${dedup[0]}"
fi

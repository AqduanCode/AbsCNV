#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

BASE=${ABSCNV_BASE:-/home/data2/DAQ/HCC1395_like_CN_50events_benchmark_20260921}
OUT="$BASE/Tools"
BAM_DIR="$BASE/alignment_cnvnator_rd_1kb_20260921/bam"
REF="$BASE/alignment_cnvnator_rd_1kb_20260921/reference/GRCh_chr1_22.fa"
CONTROL_ID=Sample_A
CASE_ID=Sample_B
CONTROL_BAM="$BAM_DIR/${CONTROL_ID}.sorted.bam"
CASE_BAM="$BAM_DIR/${CASE_ID}.sorted.bam"
SAMTOOLS=/usr/local/bin/samtools
BIN_SIZE=10000
THREADS=8
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)

CNVKIT=/home/anqiduan/miniconda3/envs/CNVkit/bin/cnvkit.py
FREEC=/home/data3/home/jjy/u251_data/G6CleanData/4-G6/controfreec/FREEC-11.6/src/freec
BICENV=/home/data3/home/jjy/u251_data/bicseq2_env
ICHOR=/home/data3/home/jjy/u251_data/ichorCNA
R_LIBS_DIR=/home/data3/home/jjy/u251_data/R_libs
ICHOR_RSCRIPT=/home/anqiduan/miniconda3/envs/ichorCNA/bin/Rscript

mkdir -p "$OUT" "$OUT/logs"
exec > >(tee -a "$OUT/logs/run_all_tools.log") 2>&1

log() { printf '[%s] %s\n' "$(date '+%F %T %z')" "$*"; }
fail() { log "ERROR: $*"; exit 1; }

for path in "$CONTROL_BAM" "$CONTROL_BAM.bai" "$CASE_BAM" "$CASE_BAM.bai" "$REF" "$REF.fai" \
  "$CNVKIT" "$FREEC" "$BICENV/bin/NBICseq-norm.pl" "$BICENV/bin/NBICseq-seg.pl" \
  "$ICHOR/scripts/runIchorCNA.R"; do
  [[ -s "$path" ]] || fail "missing or empty: $path"
done
"$SAMTOOLS" quickcheck -v "$CONTROL_BAM" "$CASE_BAM"

run_cnvkit() {
  local dir="$OUT/CNVkit_10kb"
  mkdir -p "$dir"
  log "CNVkit start"
  if [[ ! -s "$dir/${CASE_ID}.sorted.call.cns" ]]; then
    [[ -s "$dir/hg38.access.bed" ]] || "$CNVKIT" access "$REF" -o "$dir/hg38.access.bed"
    "$CNVKIT" batch "$CASE_BAM" \
      --normal "$CONTROL_BAM" --method wgs --segment-method haar \
      --target-avg-size "$BIN_SIZE" --fasta "$REF" --access "$dir/hg38.access.bed" \
      --output-reference "$dir/${CONTROL_ID}.normal_reference.10kb.cnn" \
      --output-dir "$dir" --scatter --diagram --processes "$THREADS"
  fi
  local call="$dir/${CASE_ID}.sorted.call.cns"
  [[ -s "$call" ]] || fail "CNVkit call file not produced"
  cut -f1,2,3,5,8,11 "$call" > "$dir/${CASE_ID}_vs_${CONTROL_ID}.call.selected_columns.tsv"
  cut -f1,2,3,8 "$call" > "$dir/${CASE_ID}_vs_${CONTROL_ID}.copy_number_4cols.tsv"
  log "CNVkit complete"
}

run_freec() {
  local dir="$OUT/ControlFREEC_10kb"
  mkdir -p "$dir"
  local conf="$dir/config.ini"
  cat > "$conf" <<EOF
[general]
chrLenFile = $REF.fai
ploidy = 3
breakPointThreshold = 0.8
window = $BIN_SIZE
chrFiles = $BASE/alignment_cnvnator_rd_1kb_20260921/reference/chromosomes
outputDir = $dir
numberOfProcesses = $THREADS
BedGraphOutput = FALSE
breakPointType = 1

[sample]
mateFile = $CASE_BAM
inputFormat = bam
mateOrientation = 0

[control]
mateFile = $CONTROL_BAM
inputFormat = bam
mateOrientation = 0
EOF
  log "Control-FREEC start"
  if [[ ! -s "$dir/${CASE_ID}.sorted.bam_CNVs" ]]; then
    "$FREEC" -conf "$conf" > "$dir/freec.run.log" 2>&1
  fi
  [[ -s "$dir/${CASE_ID}.sorted.bam_CNVs" ]] || fail "Control-FREEC CNV file not produced"
  log "Control-FREEC complete"
}

run_ichor() {
  local dir="$OUT/ichorCNA_10kb"
  local run_id="${CASE_ID}_vs_${CONTROL_ID}_10kb"
  local sizes="$dir/hg38_chr1-22.chrom.sizes"
  local gcwig="$dir/gc_hg38_10kb.wig"
  local control_wig="$dir/${CONTROL_ID}.reads.10kb.wig"
  local case_wig="$dir/${CASE_ID}.reads.10kb.wig"
  mkdir -p "$dir"
  "$SAMTOOLS" idxstats "$CONTROL_BAM" | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print $1"\t"$2}' > "$sizes"
  if [[ ! -s "$gcwig" ]]; then
    python3 - "$REF" "$sizes" "$gcwig.tmp.$$" "$BIN_SIZE" <<'PY'
import sys
fasta, sizes_path, output, step = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
order = [x.split('\t')[0] for x in open(sizes_path) if x.strip()]
keep, seqs, current = set(order), {c: [] for c in order}, None
for line in open(fasta):
    if line.startswith('>'):
        name = line[1:].split()[0]
        current = name if name in keep else None
    elif current is not None:
        seqs[current].append(line.strip().upper())
with open(output, 'w') as out:
    for chrom in order:
        seq = ''.join(seqs[chrom])
        out.write(f'fixedStep chrom={chrom} start=1 step={step} span={step}\n')
        for start in range(0, len(seq), step):
            window = seq[start:start+step]
            valid = sum(base in 'ACGT' for base in window)
            out.write('-1\n' if valid == 0 else f'{sum(base in "GC" for base in window)/valid:.6f}\n')
PY
    mv "$gcwig.tmp.$$" "$gcwig"
  fi
  make_wig() {
    local bam="$1" wig="$2"
    [[ -s "$wig" ]] && return 0
    "$SAMTOOLS" view -F 3844 "$bam" "${CHROMS[@]}" | awk -v BS="$BIN_SIZE" -v SIZES="$sizes" '
      BEGIN {while ((getline line < SIZES)>0){split(line,f,"\t"); n++; chr[n]=f[1]; bins[f[1]]=int((f[2]+BS-1)/BS)}}
      {b=int(($4-1)/BS)+1; counts[$3 SUBSEP b]++}
      END {for(i=1;i<=n;i++){c=chr[i]; print "fixedStep chrom="c" start=1 step="BS" span="BS; for(b=1;b<=bins[c];b++) print counts[c SUBSEP b]+0}}' > "$wig.tmp.$$"
    mv "$wig.tmp.$$" "$wig"
  }
  log "ichorCNA read counting start"
  make_wig "$CONTROL_BAM" "$control_wig" & local p1=$!
  make_wig "$CASE_BAM" "$case_wig" & local p2=$!
  wait "$p1"; wait "$p2"
  log "ichorCNA model start"
  if [[ ! -s "$dir/${run_id}.seg" ]]; then
    R_LIBS="$R_LIBS_DIR" "$ICHOR_RSCRIPT" "$ICHOR/scripts/runIchorCNA.R" \
      --id "$run_id" --WIG "$case_wig" --NORMWIG "$control_wig" \
      --gcWig "$gcwig" --mapWig None --normalPanel None \
      --ploidy 'c(2,3,4)' --normal 'c(0,0.1,0.2,0.3,0.4,0.5)' --maxCN 7 \
      --includeHOMD FALSE --chrs 'c(1:22)' --chrTrain 'c(1:22)' --chrNormalize 'c(1:22)' \
      --estimateNormal TRUE --estimatePloidy TRUE --estimateScPrevalence TRUE --scStates 'c(1,3)' \
      --centromere "$ICHOR/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt" \
      --genomeBuild hg38 --genomeStyle UCSC --exons.bed None \
      --txnE 0.9999999 --txnStrength 10000000 --plotFileType pdf --plotYLim 'c(-2,4)' \
      --outDir "$dir" > "$dir/${run_id}.log" 2>&1
  fi
  [[ -s "$dir/${run_id}.seg" ]] || fail "ichorCNA segment file not produced"
  log "ichorCNA complete"
}

run_bicseq2() {
  local dir="$OUT/BICseq2_10kb"
  local read_len=100 frag_size=200 lambda=2
  export PATH="$BICENV/bin:$PATH"
  mkdir -p "$dir"/{chrom_fasta,mappability_nonN,read_positions/{Sample_A,Sample_B},normalized_bins_10kb/{Sample_A,Sample_B},tmp_10kb/{Sample_A,Sample_B,seg},results}
  log "BICseq2 reference preparation start"
  for chrom in "${CHROMS[@]}"; do
    local fa="$dir/chrom_fasta/$chrom.fa" map="$dir/mappability_nonN/$chrom.nonN.map"
    [[ -s "$fa" ]] || "$SAMTOOLS" faidx "$REF" "$chrom" > "$fa"
    if [[ ! -s "$map" ]]; then
      perl -ne 'BEGIN{$pos=0;$start=0} next if /^>/; chomp; for $b(split //){$pos++; if($b=~/[Nn]/){if($start){print $start,"\t",$pos-1,"\n";$start=0}}else{$start||=$pos}} END{print $start,"\t",$pos,"\n" if $start}' "$fa" > "$map"
    fi
  done
  extract_positions() {
    local id="$1" bam="$2"
    for chrom in "${CHROMS[@]}"; do
      local seq="$dir/read_positions/$id/$chrom.seq"
      [[ -s "$seq" ]] || "$SAMTOOLS" view -f 64 -F 260 "$bam" "$chrom" | awk '{print $4}' > "$seq"
    done
  }
  run_norm() {
    local id="$1"
    local conf="$dir/config_norm_${id}_10kb.tsv"
    local params="$dir/results/${id}.10kb.norm_params.txt"
    { echo -e 'chrom\tfa\tmappability\tread_positions\tnormalized_bins'; for chrom in "${CHROMS[@]}"; do echo -e "$chrom\t$dir/chrom_fasta/$chrom.fa\t$dir/mappability_nonN/$chrom.nonN.map\t$dir/read_positions/$id/$chrom.seq\t$dir/normalized_bins_10kb/$id/$chrom.norm.bin"; done; } > "$conf"
    NBICseq-norm.pl -l="$read_len" -s="$frag_size" -b="$BIN_SIZE" --tmp="$dir/tmp_10kb/$id" "$conf" "$params"
  }
  log "BICseq2 read position extraction start"
  extract_positions "$CONTROL_ID" "$CONTROL_BAM"
  extract_positions "$CASE_ID" "$CASE_BAM"
  log "BICseq2 normalization start"
  run_norm "$CONTROL_ID"; run_norm "$CASE_ID"
  local conf="$dir/config_seg_${CASE_ID}_vs_${CONTROL_ID}_10kb.tsv"
  local result="$dir/results/${CASE_ID}_vs_${CONTROL_ID}.10kb.bicseq2_segments.tsv"
  { echo -e 'chrom\tcase_normalized_bins\tcontrol_normalized_bins'; for chrom in "${CHROMS[@]}"; do echo -e "$chrom\t$dir/normalized_bins_10kb/$CASE_ID/$chrom.norm.bin\t$dir/normalized_bins_10kb/$CONTROL_ID/$chrom.norm.bin"; done; } > "$conf"
  log "BICseq2 segmentation start"
  NBICseq-seg.pl --control --lambda="$lambda" --tmp="$dir/tmp_10kb/seg" \
    --fig="$dir/results/${CASE_ID}_vs_${CONTROL_ID}.10kb.profile.png" \
    --title="$CASE_ID vs $CONTROL_ID BICseq2 CNV 10kb" "$conf" "$result"
  [[ -s "$result" ]] || fail "BICseq2 segment file not produced"
  awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++)if($i=="log2.copyRatio")c=i; print $0,"copyRatio","status"; next}{r=2^($c+0); s=(r<0.9?"loss":(r>1.1?"gain":"neutral")); print $0,r,s}' "$result" > "$dir/results/${CASE_ID}_vs_${CONTROL_ID}.10kb.bicseq2_segments.annotated.tsv"
  log "BICseq2 complete"
}

case "${RUN_ONLY:-all}" in
  ichor)
    log "ichorCNA-only workflow started"
    run_ichor
    log "ichorCNA-only workflow completed"
    ;;
  bicseq2)
    log "BICseq2-only workflow started"
    run_bicseq2
    log "BICseq2-only workflow completed"
    ;;
  all)
    log "Four-tool workflow started: case=$CASE_ID control=$CONTROL_ID bin=$BIN_SIZE"
    run_cnvkit
    run_freec
    run_ichor
    run_bicseq2
    date '+%F %T %z' > "$OUT/COMPLETED.txt"
    log "Four-tool workflow completed"
    ;;
  *)
    fail "RUN_ONLY must be one of: all, ichor, bicseq2"
    ;;
esac

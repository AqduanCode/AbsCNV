#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

PROJECT=${ABSCNV_PROJECT:-/home/data2/DAQ/HCC1395_like_Sample_A_single_sample_CNV_benchmark_20261004}
OUT="$PROJECT/Tools/ControlFREEC_10kb"
BAM="$PROJECT/alignment/bam/Sample_A.sorted.bam"
REF="$PROJECT/reference/GRCh_chr1_22.fa"
FREEC=/home/data3/home/jjy/u251_data/G6CleanData/4-G6/controfreec/FREEC-11.6/src/freec
SAMTOOLS=/usr/local/bin/samtools
STATUS="$PROJECT/PIPELINE_STATUS.tsv"
ALIGN_PID_FILE="$PROJECT/logs/01_align_sampleA.pid"

mkdir -p "$OUT"
record() {
  printf '%s\t%s\t%s\t%s\n' "$(date '+%F %T %z')" "$1" "$2" "$3" >> "$STATUS"
}

while [[ ! -e "$PROJECT/alignment/qc/ALIGNMENT_COMPLETE" ]]; do
  if [[ -s "$ALIGN_PID_FILE" ]] && ! kill -0 "$(cat "$ALIGN_PID_FILE")" 2>/dev/null; then
    record freec failed "alignment stopped before completion"
    exit 1
  fi
  sleep 60
done

[[ -x "$FREEC" ]] || { record freec failed "missing executable: $FREEC"; exit 1; }
[[ -s "$REF.fai" ]] || { record freec failed "missing reference index: $REF.fai"; exit 1; }
[[ -s "$BAM" && -s "$BAM.bai" ]] || { record freec failed "missing BAM or index"; exit 1; }
"$SAMTOOLS" quickcheck -v "$BAM"

cat > "$OUT/config.ini" <<EOF
[general]
chrLenFile = $REF.fai
ploidy = 2
breakPointThreshold = 0.8
breakPointType = 4
window = 10000
chrFiles = $PROJECT/reference/chromosomes
outputDir = $OUT
numberOfProcesses = 16
BedGraphOutput = FALSE

[sample]
mateFile = $BAM
inputFormat = BAM
mateOrientation = 0
EOF

rm -f "$OUT/COMPLETED.txt"
record freec started "single-sample WGS; 10 kb; ploidy 2; overwrite previous ploidy-3 outputs"
"$FREEC" -conf "$OUT/config.ini" > "$OUT/freec.run.log" 2>&1

[[ -s "$OUT/Sample_A.sorted.bam_CNVs" ]] || {
  record freec failed "CNV file missing; see $OUT/freec.run.log"
  exit 1
}
record freec complete "calls=$OUT/Sample_A.sorted.bam_CNVs"
touch "$OUT/COMPLETED.txt"

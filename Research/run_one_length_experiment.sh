#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

LABEL=$1
PROJECT=${ABSCNV_PROJECT:-/home/data2/DAQ/HCC1395_like_CN_length_separate_benchmark_20260928}
EXP="$PROJECT/experiments/$LABEL"
mkdir -p "$EXP/logs" "$EXP/qc"

printf '[%s] %s simulation starting\n' "$(date '+%F %T %z')" "$LABEL"
bash "$RESEARCH_DIR/generate_length_benchmark_fastq.sh" "$LABEL"
printf '[%s] %s waiting for alignment lock\n' "$(date '+%F %T %z')" "$LABEL"

# Six simulations run concurrently, but alignment is serialized to limit
# temporary BAM space and avoid six simultaneous bwa/samtools I/O streams.
flock "$PROJECT/qc/alignment.lock" \
  bash "$RESEARCH_DIR/align_length_benchmark_reads.sh" "$LABEL"

date '+%F %T %z' > "$EXP/qc/full_pipeline_completed_at.txt"
printf '[%s] %s full pipeline complete\n' "$(date '+%F %T %z')" "$LABEL"

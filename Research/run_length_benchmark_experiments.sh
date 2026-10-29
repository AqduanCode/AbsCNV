#!/usr/bin/env bash
set -euo pipefail
export RESEARCH_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

PROJECT=${ABSCNV_PROJECT:-/home/data2/DAQ/HCC1395_like_CN_length_separate_benchmark_20260928}
REFERENCE_PROJECT=${ABSCNV_LENGTH_SOURCE:-$PROJECT/simulation_sources}
REF_DIR=/home/data3/home/jjy/CNVsim/refs
mkdir -p "$PROJECT/logs" "$PROJECT/qc" "$PROJECT/experiments"

available_bytes=$(df --output=avail -B1 "$PROJECT" | tail -n 1 | tr -d ' ')
if (( available_bytes < 650000000000 )); then
  echo "ERROR: less than 650 GB free before six-length pipeline: $available_bytes" >&2
  exit 1
fi

python3 "$RESEARCH_DIR/build_length_benchmark_profiles.py" \
  "$REFERENCE_PROJECT/cnv_gain_cn_median.txt" \
  "$REFERENCE_PROJECT/cnv_loss_cn_median.txt" \
  "$REF_DIR" \
  "$REFERENCE_PROJECT/Sample_B_vs_Sample_A.master_120_events.tsv" \
  "$PROJECT" \
  "$REFERENCE_PROJECT/Sample_A.absolute_CN_segments.tsv"

printf 'length_group\tpid\tlog\n' > "$PROJECT/qc/controller_jobs.tsv"
for label in 20kb 50kb 100kb 200kb 500kb 1000kb; do
  mkdir -p "$PROJECT/experiments/$label/logs"
  nohup bash "$RESEARCH_DIR/run_one_length_experiment.sh" "$label" \
    > "$PROJECT/experiments/$label/logs/full_pipeline.nohup.log" 2>&1 < /dev/null &
  pid=$!
  printf '%s\t%s\t%s\n' "$label" "$pid" "$PROJECT/experiments/$label/logs/full_pipeline.nohup.log" \
    >> "$PROJECT/qc/controller_jobs.tsv"
done

date '+%F %T %z' > "$PROJECT/qc/launched_at.txt"
echo "Launched six independent length workflows"
cat "$PROJECT/qc/controller_jobs.tsv"

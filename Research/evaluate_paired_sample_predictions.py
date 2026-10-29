#!/usr/bin/env python3
from __future__ import annotations

import csv
import argparse
import importlib.util
import math
from pathlib import Path


REFERENCE_SCRIPT = Path(__file__).with_name("paired_interval_metrics.py")
TRUTH_PATH = NEW_ABSCNV_PATH = OTHER_CALLS_DIR = OUTPUT_DIR = None

EVENT_STANDARDS = {
    "single80_one_to_one": (0.80, "either", "one_to_one"),
    "reciprocal80_one_to_one": (0.80, "reciprocal", "one_to_one"),
    "single80_many_to_many": (0.80, "either", "many_to_many"),
}


def load_reference_module():
    spec = importlib.util.spec_from_file_location("cnv_evaluation_reference", REFERENCE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load reference evaluation script: {REFERENCE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def metric(matched_predictions, predicted_events, matched_truth, truth_events):
    precision = matched_predictions / predicted_events if predicted_events else math.nan
    recall = matched_truth / truth_events if truth_events else math.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def many_to_many_eval(evaluation, tool, predictions, truths, standard):
    candidates = evaluation.candidates(predictions, truths, 0.80, "either")
    matches = []
    for _, iou, overlap_bp, prediction_fraction, truth_fraction, prediction, truth in candidates:
        matches.append(
            dict(
                standard=standard,
                tool=tool,
                prediction_id=prediction["id"],
                truth_id=truth["id"],
                chrom=prediction["chrom"],
                direction=prediction["direction"],
                prediction_start=prediction["start"],
                prediction_end=prediction["end"],
                truth_start=truth["start"],
                truth_end=truth["end"],
                overlap_bp=overlap_bp,
                prediction_overlap_fraction=prediction_fraction,
                truth_overlap_fraction=truth_fraction,
                iou=iou,
            )
        )

    rows = []
    for direction in ("all", "gain", "loss"):
        selected_predictions = [p for p in predictions if direction == "all" or p["direction"] == direction]
        selected_truths = [t for t in truths if direction == "all" or t["direction"] == direction]
        selected_matches = [m for m in matches if direction == "all" or m["direction"] == direction]
        matched_prediction_ids = {m["prediction_id"] for m in selected_matches}
        matched_truth_ids = {m["truth_id"] for m in selected_matches}
        precision, recall, f1 = metric(
            len(matched_prediction_ids),
            len(selected_predictions),
            len(matched_truth_ids),
            len(selected_truths),
        )
        rows.append(
            dict(
                standard=standard,
                tool=tool,
                direction=direction,
                matched_predictions=len(matched_prediction_ids),
                matched_truth_events=len(matched_truth_ids),
                false_positive_predictions=len(selected_predictions) - len(matched_prediction_ids),
                false_negative_truth_events=len(selected_truths) - len(matched_truth_ids),
                predicted_events=len(selected_predictions),
                truth_events=len(selected_truths),
                precision=precision,
                recall=recall,
                f1=f1,
            )
        )

    used_predictions = {m["prediction_id"] for m in matches}
    used_truths = {m["truth_id"] for m in matches}
    unmatched_predictions = [dict(standard=standard, tool=tool, **p) for p in predictions if p["id"] not in used_predictions]
    unmatched_truth = [dict(standard=standard, tool=tool, **t) for t in truths if t["id"] not in used_truths]
    return rows, matches, unmatched_predictions, unmatched_truth


def main():
    global TRUTH_PATH, NEW_ABSCNV_PATH, OTHER_CALLS_DIR, OUTPUT_DIR
    parser = argparse.ArgumentParser(description="Original four evaluation standards; configurable archived inputs")
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    TRUTH_PATH = args.dataset / "benchmark/Sample_B_vs_Sample_A.benchmark.tsv"
    NEW_ABSCNV_PATH = args.dataset / "output/CNVresult_filter.csv"
    OTHER_CALLS_DIR = args.dataset / "normalized_calls"
    OUTPUT_DIR = args.output_dir
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    evaluation = load_reference_module()
    truths = evaluation.load_truth(TRUTH_PATH)
    tools = {
        "AbsCNV_background_updated": evaluation.load_my(NEW_ABSCNV_PATH),
    }
    for path in sorted(OTHER_CALLS_DIR.glob("*.tsv")):
        tools[path.stem] = evaluation.load_normalized(path)

    evaluation.validate(truths, "truth")
    for tool, predictions in tools.items():
        evaluation.validate(predictions, tool)

    event_rows = []
    matches = []
    unmatched_predictions = []
    unmatched_truth = []
    summary = []

    for standard, (threshold, mode, multiplicity) in EVENT_STANDARDS.items():
        for tool, predictions in tools.items():
            if multiplicity == "many_to_many":
                rows, selected, unpred, untruth = many_to_many_eval(
                    evaluation, tool, predictions, truths, standard
                )
            else:
                rows, selected, unpred, untruth = evaluation.event_eval(
                    tool, predictions, truths, standard, threshold, mode, multiplicity
                )
            event_rows.extend(rows)
            matches.extend(selected)
            unmatched_predictions.extend(unpred)
            unmatched_truth.extend(untruth)
            all_row = next(row for row in rows if row["direction"] == "all")
            summary.append(
                dict(
                    criterion=standard,
                    metric_level="event",
                    tool=tool,
                    precision=all_row["precision"],
                    recall=all_row["recall"],
                    f1=all_row["f1"],
                    predicted_events=all_row["predicted_events"],
                    truth_events=all_row["truth_events"],
                    matched_predictions=all_row["matched_predictions"],
                    matched_truth_events=all_row["matched_truth_events"],
                )
            )

    base_rows = []
    for tool, predictions in tools.items():
        rows = evaluation.base_eval(tool, predictions, truths, "base_pair")
        base_rows.extend(rows)
        all_row = next(row for row in rows if row["direction"] == "all")
        summary.append(
            dict(
                criterion="base_pair",
                metric_level="base",
                tool=tool,
                precision=all_row["precision"],
                recall=all_row["recall"],
                f1=all_row["f1"],
                predicted_events="",
                truth_events="",
                matched_predictions="",
                matched_truth_events="",
            )
        )

    evaluation.write_table(OUTPUT_DIR / "four_standard_f1_summary.csv", summary)
    evaluation.write_table(OUTPUT_DIR / "four_standard_event_metrics.csv", event_rows)
    evaluation.write_table(OUTPUT_DIR / "four_standard_base_metrics.csv", base_rows)
    evaluation.write_table(OUTPUT_DIR / "four_standard_event_matches.csv", matches)
    evaluation.write_table(OUTPUT_DIR / "four_standard_unmatched_predictions.csv", unmatched_predictions)
    evaluation.write_table(OUTPUT_DIR / "four_standard_unmatched_truth.csv", unmatched_truth)

    definition = f"""Four F1 evaluation criteria

Truth: {TRUTH_PATH}
Updated AbsCNV: {NEW_ABSCNV_PATH}
Coordinates: 0-based half-open; length=end-start.
Direction: Sample_B relative to Sample_A; gain/loss must agree.
1. Single-direction 80%, one-to-one: prediction OR truth overlap fraction >=80%; greedy best one-to-one matching.
2. Reciprocal 80%, one-to-one: prediction AND truth overlap fractions >=80%; greedy best one-to-one matching.
3. Single-direction 80%, many-to-many: every direction-consistent pair satisfying prediction OR truth overlap fraction >=80% is accepted. Precision counts unique matched predictions; recall counts unique matched truths.
4. Base-pair level: direction-specific interval unions and same-direction intersections.
"""
    (OUTPUT_DIR / "four_standard_evaluation_definition.txt").write_text(definition, encoding="utf-8")

    for row in summary:
        print(f"{row['criterion']}\t{row['tool']}\t{row['f1']:.6f}")


if __name__ == "__main__":
    main()

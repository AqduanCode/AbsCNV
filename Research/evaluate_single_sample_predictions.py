#!/usr/bin/env python3
"""Evaluate the latest five single-sample CNV outputs by paired-sample rules.

The four criteria reproduce the definitions and greedy event matching used in
the paired-sample semi-simulation evaluation under 背景更新后的结果/半模拟数据/evaluate.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


EVENT_STANDARDS = {
    "single80_one_to_one": (0.80, "either", "one_to_one"),
    "reciprocal80_one_to_one": (0.80, "reciprocal", "one_to_one"),
    "single80_many_to_many": (0.80, "either", "many_to_many"),
}
TOOLS = ("MyMethod", "CNVnator", "ControlFREEC", "RDXplorer", "CONY")


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def record(call) -> dict:
    return {
        "id": call.call_id,
        "chrom": call.chrom,
        "start": call.start,
        "end": call.end,
        "direction": call.direction,
        "source": call.source,
    }


def validate(records: list[dict], label: str) -> None:
    identifiers = set()
    for row in records:
        if row["id"] in identifiers:
            raise ValueError(f"{label}: duplicate id {row['id']}")
        identifiers.add(row["id"])
        if row["direction"] not in {"gain", "loss"}:
            raise ValueError(f"{label}: invalid direction {row['direction']}")
        if row["start"] >= row["end"]:
            raise ValueError(f"{label}: non-positive interval {row}")


def overlap(left: dict, right: dict) -> int:
    return max(0, min(left["end"], right["end"]) - max(left["start"], right["start"]))


def metric(tp: int, predicted: int, truth: int):
    precision = tp / predicted if predicted else math.nan
    recall = tp / truth if truth else math.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def candidates(predictions: list[dict], truths: list[dict], threshold: float,
               mode: str):
    truth_groups = defaultdict(list)
    for truth in truths:
        truth_groups[(truth["chrom"], truth["direction"])].append(truth)
    answer = []
    for prediction in predictions:
        prediction_length = prediction["end"] - prediction["start"]
        for truth in truth_groups[(prediction["chrom"], prediction["direction"])]:
            overlap_bp = overlap(prediction, truth)
            if overlap_bp <= 0:
                continue
            truth_length = truth["end"] - truth["start"]
            prediction_fraction = overlap_bp / prediction_length
            truth_fraction = overlap_bp / truth_length
            accepted = (
                prediction_fraction >= threshold and truth_fraction >= threshold
                if mode == "reciprocal"
                else prediction_fraction >= threshold or truth_fraction >= threshold
            )
            if accepted:
                iou = overlap_bp / (prediction_length + truth_length - overlap_bp)
                answer.append((
                    min(prediction_fraction, truth_fraction), iou, overlap_bp,
                    prediction_fraction, truth_fraction, prediction, truth,
                ))
    return sorted(answer, key=lambda item: (item[0], item[1], item[2]), reverse=True)


def event_eval(tool: str, predictions: list[dict], truths: list[dict], standard: str,
               threshold: float, mode: str, multiplicity: str):
    candidate_pairs = candidates(predictions, truths, threshold, mode)
    used_predictions = set()
    used_truths = set()
    matches = []

    if multiplicity == "many_to_many":
        # Preserve every qualifying pair. Precision and recall later count
        # unique prediction and truth identifiers separately.
        selected = candidate_pairs
    else:
        selected = []
        for item in candidate_pairs:
            prediction, truth = item[-2], item[-1]
            if prediction["id"] in used_predictions or truth["id"] in used_truths:
                continue
            used_predictions.add(prediction["id"])
            used_truths.add(truth["id"])
            selected.append(item)

    for score, iou, overlap_bp, prediction_fraction, truth_fraction, prediction, truth in selected:
        used_predictions.add(prediction["id"])
        used_truths.add(truth["id"])
        matches.append({
            "standard": standard,
            "tool": tool,
            "prediction_id": prediction["id"],
            "truth_id": truth["id"],
            "chrom": prediction["chrom"],
            "direction": prediction["direction"],
            "prediction_start": prediction["start"],
            "prediction_end": prediction["end"],
            "truth_start": truth["start"],
            "truth_end": truth["end"],
            "overlap_bp": overlap_bp,
            "prediction_overlap_fraction": prediction_fraction,
            "truth_overlap_fraction": truth_fraction,
            "iou": iou,
        })

    rows = []
    for direction in ("all", "gain", "loss"):
        selected_predictions = [
            row for row in predictions if direction == "all" or row["direction"] == direction
        ]
        selected_truths = [
            row for row in truths if direction == "all" or row["direction"] == direction
        ]
        selected_matches = [
            row for row in matches if direction == "all" or row["direction"] == direction
        ]
        matched_prediction_ids = {row["prediction_id"] for row in selected_matches}
        matched_truth_ids = {row["truth_id"] for row in selected_matches}
        if multiplicity == "one_to_one":
            precision, recall, f1 = metric(
                len(matched_prediction_ids), len(selected_predictions), len(selected_truths)
            )
        else:
            precision = (len(matched_prediction_ids) / len(selected_predictions)
                         if selected_predictions else math.nan)
            recall = (len(matched_truth_ids) / len(selected_truths)
                      if selected_truths else math.nan)
            f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append({
            "standard": standard,
            "tool": tool,
            "direction": direction,
            "matched_predictions": len(matched_prediction_ids),
            "matched_truth_events": len(matched_truth_ids),
            "false_positive_predictions": len(selected_predictions) - len(matched_prediction_ids),
            "false_negative_truth_events": len(selected_truths) - len(matched_truth_ids),
            "predicted_events": len(selected_predictions),
            "truth_events": len(selected_truths),
            "precision": precision,
            "recall": recall,
            "f1": f1,
        })

    unmatched_predictions = [
        {"standard": standard, "tool": tool, **row}
        for row in predictions if row["id"] not in used_predictions
    ]
    unmatched_truths = [
        {"standard": standard, "tool": tool, **row}
        for row in truths if row["id"] not in used_truths
    ]
    return rows, matches, unmatched_predictions, unmatched_truths


def merge_by_key(records: list[dict]):
    groups = defaultdict(list)
    for row in records:
        groups[(row["chrom"], row["direction"])].append((row["start"], row["end"]))
    merged_groups = {}
    for key, intervals in groups.items():
        merged = []
        for start, end in sorted(intervals):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        merged_groups[key] = merged
    return merged_groups


def sum_bp(groups, direction=None) -> int:
    return sum(
        end - start
        for (_, row_direction), intervals in groups.items()
        if direction is None or row_direction == direction
        for start, end in intervals
    )


def intersect_bp(left, right, direction=None) -> int:
    total = 0
    for key, left_intervals in left.items():
        if direction is not None and key[1] != direction:
            continue
        right_intervals = right.get(key, [])
        i = j = 0
        while i < len(left_intervals) and j < len(right_intervals):
            total += max(
                0,
                min(left_intervals[i][1], right_intervals[j][1])
                - max(left_intervals[i][0], right_intervals[j][0]),
            )
            if left_intervals[i][1] <= right_intervals[j][1]:
                i += 1
            else:
                j += 1
    return total


def base_eval(tool: str, predictions: list[dict], truths: list[dict]):
    prediction_map = merge_by_key(predictions)
    truth_map = merge_by_key(truths)
    rows = []
    for label, direction in (("all", None), ("gain", "gain"), ("loss", "loss")):
        predicted_bp = sum_bp(prediction_map, direction)
        truth_bp = sum_bp(truth_map, direction)
        matched_bp = intersect_bp(prediction_map, truth_map, direction)
        precision, recall, f1 = metric(matched_bp, predicted_bp, truth_bp)
        rows.append({
            "standard": "base_pair",
            "tool": tool,
            "direction": label,
            "matched_bp": matched_bp,
            "predicted_cnv_bp": predicted_bp,
            "truth_cnv_bp": truth_bp,
            "false_positive_bp": predicted_bp - matched_bp,
            "false_negative_bp": truth_bp - matched_bp,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base", type=Path,
        required=True,
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    base = args.base
    output_dir = args.output_dir

    # Reuse the already-validated latest raw parsers without changing calls.
    from read_single_sample_predictions import (  # noqa: E402
        read_cnvnator,
        read_cony,
        read_freec,
        read_mymethod,
        read_rdxplorer,
    )
    from read_single_sample_predictions import read_event_tsv  # noqa: E402

    other = base / "othertools"
    input_paths = {
        "truth": base / "benchmark" / "Sample_A.truth_events.tsv",
        "MyMethod": other / "MyMethod_raw.tsv",
        "CNVnator": other / "CNVnator_raw.tsv",
        "ControlFREEC": other / "ControlFREEC_raw.tsv",
        "RDXplorer": other / "RDXplorer_raw.tsv",
        "CONY": other / "CONY_raw.tsv",
    }
    missing = [str(path) for path in input_paths.values()
               if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing or empty inputs:\n" + "\n".join(missing))

    truths = [record(call) for call in read_event_tsv(input_paths["truth"])]
    calls_by_tool = {
        "MyMethod": [record(call) for call in read_mymethod(input_paths["MyMethod"])],
        "CNVnator": [record(call) for call in read_cnvnator(input_paths["CNVnator"])],
        "ControlFREEC": [record(call) for call in read_freec(input_paths["ControlFREEC"])],
        "RDXplorer": [record(call) for call in read_rdxplorer(input_paths["RDXplorer"])],
        "CONY": [record(call) for call in read_cony(input_paths["CONY"])],
    }
    for tool, calls in calls_by_tool.items():
        write_csv(output_dir / "standardized_calls" / (tool + ".csv"), calls,
                  ["id", "chrom", "start", "end", "direction", "source"])
    validate(truths, "truth")
    for tool, calls in calls_by_tool.items():
        validate(calls, tool)

    event_rows = []
    match_rows = []
    unmatched_prediction_rows = []
    unmatched_truth_rows = []
    summary_rows = []

    for standard, (threshold, mode, multiplicity) in EVENT_STANDARDS.items():
        for tool in TOOLS:
            rows, matches, unmatched_predictions, unmatched_truths = event_eval(
                tool, calls_by_tool[tool], truths, standard, threshold, mode, multiplicity
            )
            event_rows.extend(rows)
            match_rows.extend(matches)
            unmatched_prediction_rows.extend(unmatched_predictions)
            unmatched_truth_rows.extend(unmatched_truths)
            all_row = next(row for row in rows if row["direction"] == "all")
            summary_rows.append({
                "criterion": standard,
                "metric_level": "event",
                "tool": tool,
                "precision": all_row["precision"],
                "recall": all_row["recall"],
                "f1": all_row["f1"],
                "predicted_events": all_row["predicted_events"],
                "truth_events": all_row["truth_events"],
                "matched_predictions": all_row["matched_predictions"],
                "matched_truth_events": all_row["matched_truth_events"],
            })

    base_rows = []
    for tool in TOOLS:
        rows = base_eval(tool, calls_by_tool[tool], truths)
        base_rows.extend(rows)
        all_row = next(row for row in rows if row["direction"] == "all")
        summary_rows.append({
            "criterion": "base_pair",
            "metric_level": "base",
            "tool": tool,
            "precision": all_row["precision"],
            "recall": all_row["recall"],
            "f1": all_row["f1"],
            "predicted_events": "",
            "truth_events": "",
            "matched_predictions": "",
            "matched_truth_events": "",
        })

    write_csv(output_dir / "single_sample_four_standard_summary.csv", summary_rows, [
        "criterion", "metric_level", "tool", "precision", "recall", "f1",
        "predicted_events", "truth_events", "matched_predictions", "matched_truth_events",
    ])
    write_csv(output_dir / "single_sample_event_metrics.csv", event_rows, [
        "standard", "tool", "direction", "matched_predictions", "matched_truth_events",
        "false_positive_predictions", "false_negative_truth_events", "predicted_events",
        "truth_events", "precision", "recall", "f1",
    ])
    write_csv(output_dir / "single_sample_base_metrics.csv", base_rows, [
        "standard", "tool", "direction", "matched_bp", "predicted_cnv_bp",
        "truth_cnv_bp", "false_positive_bp", "false_negative_bp", "precision",
        "recall", "f1",
    ])
    write_csv(output_dir / "single_sample_event_matches.csv", match_rows, [
        "standard", "tool", "prediction_id", "truth_id", "chrom", "direction",
        "prediction_start", "prediction_end", "truth_start", "truth_end", "overlap_bp",
        "prediction_overlap_fraction", "truth_overlap_fraction", "iou",
    ])
    write_csv(output_dir / "single_sample_unmatched_predictions.csv",
              unmatched_prediction_rows,
              ["standard", "tool", "id", "chrom", "start", "end", "direction", "source"])
    write_csv(output_dir / "single_sample_unmatched_truth.csv", unmatched_truth_rows,
              ["standard", "tool", "id", "chrom", "start", "end", "direction", "source"])

    definitions = f"""Single-sample evaluation using the paired-sample semi-simulation criteria

Truth: {input_paths['truth']}
Inputs: latest top-level othertools/*_raw.tsv files
Coordinates: normalized to 0-based half-open; length=end-start.
Direction: gain/loss must agree.
1. Single-direction 80%, one-to-one: prediction OR truth overlap fraction >=80%; greedy best one-to-one matching.
2. Reciprocal 80%, one-to-one: prediction AND truth overlap fractions >=80%; greedy best one-to-one matching.
3. Single-direction 80%, many-to-many: every direction-consistent qualifying pair is accepted. Precision counts unique matched predictions; recall counts unique matched truths.
4. Base-pair level: direction-specific interval unions and same-direction intersections.
"""
    (output_dir / "single_sample_evaluation_definitions.txt").write_text(
        definitions, encoding="utf-8"
    )

    for row in summary_rows:
        print(
            f"{row['criterion']}\t{row['tool']}\t"
            f"{float(row['precision']):.9f}\t{float(row['recall']):.9f}\t"
            f"{float(row['f1']):.9f}"
        )


if __name__ == "__main__":
    main()

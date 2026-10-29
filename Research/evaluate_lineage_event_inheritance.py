#!/usr/bin/env python3
"""Plot P2 fates from overlaps between filtered P0-P1 and P1-P2 CNV events."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


SAMPLE_COLUMNS = {
    "U1": "G6",
    "U1.1": "G6A7first", "U1.2": "G6B11first", "U1.3": "G6C7first",
    "U1.1.1": "G6A7E1", "U1.1.2": "G6A7E11",
    "U1.2.1": "G6B11C1", "U1.2.2": "G6B11C6",
    "U1.3.1": "G6C7E1", "U1.3.2": "G6C7E9",
}
BRANCHES = {
    "U1.1": ("U1.1.1", "U1.1.2"),
    "U1.2": ("U1.2.1", "U1.2.2"),
    "U1.3": ("U1.3.1", "U1.3.2"),
}
CATEGORY_ORDER = (
    "Exact inheritance",
    "Further same-direction change",
    "Opposite-direction change",
)
COLORS = {
    "Exact inheritance": "#4E79A7",
    "Further same-direction change": "#F2A541",
    "Opposite-direction change": "#8E6C8A",
}
OUTPUT_STEM = "six_true_inheritance_three_fates_P1_event_overlap"


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def oriented_pair_events(cnv: pd.DataFrame, sample_a: str, sample_b: str) -> pd.DataFrame:
    """Return pair events oriented as sample_a -> sample_b."""
    source_a = SAMPLE_COLUMNS[sample_a]
    source_b = SAMPLE_COLUMNS[sample_b]
    forward = cnv[
        (cnv["pair_sampleA_name"] == source_a)
        & (cnv["pair_sampleB_name"] == source_b)
    ].copy()
    forward["cn_a"] = pd.to_numeric(forward["copy_number_A"], errors="coerce")
    forward["cn_b"] = pd.to_numeric(forward["copy_number_B"], errors="coerce")
    forward["source_orientation"] = "forward"

    reverse = cnv[
        (cnv["pair_sampleA_name"] == source_b)
        & (cnv["pair_sampleB_name"] == source_a)
    ].copy()
    reverse["cn_a"] = pd.to_numeric(reverse["copy_number_B"], errors="coerce")
    reverse["cn_b"] = pd.to_numeric(reverse["copy_number_A"], errors="coerce")
    reverse["source_orientation"] = "reverse"

    events = pd.concat([forward, reverse], ignore_index=True)
    if events.empty:
        raise ValueError(f"No filtered events found for {sample_a}-{sample_b}")
    for col in ("start_pos", "end_pos", "length"):
        events[col] = pd.to_numeric(events[col], errors="raise").astype(np.int64)
    if events[["cn_a", "cn_b"]].isna().any().any():
        raise ValueError(f"Missing copy-number state in {sample_a}-{sample_b}")
    if (events["end_pos"] <= events["start_pos"]).any():
        raise ValueError(f"Non-positive interval in {sample_a}-{sample_b}")
    calculated_length = events["end_pos"] - events["start_pos"]
    if not (calculated_length == events["length"]).all():
        raise ValueError(f"The length column is inconsistent with end_pos-start_pos in {sample_a}-{sample_b}")
    if (events["cn_a"] == events["cn_b"]).any():
        raise ValueError(f"Equal CN states found in filtered events for {sample_a}-{sample_b}")

    events = events.sort_values(["chr", "start_pos", "end_pos"]).reset_index(drop=True)
    events["event_id"] = [f"{sample_a}_{sample_b}_HC{i + 1}" for i in range(len(events))]
    for chrom, group in events.groupby("chr", sort=False):
        previous_end = None
        for row in group.itertuples(index=False):
            if previous_end is not None and int(row.start_pos) < previous_end:
                raise ValueError(f"Overlapping filtered events in {sample_a}-{sample_b} on {chrom}")
            previous_end = int(row.end_pos)
    return events


def change_category(p0_cn: float, p1_cn: float, pair_p1_cn: float, p2_cn: float) -> str:
    p0_to_p1 = p1_cn - p0_cn
    p1_to_p2 = p2_cn - pair_p1_cn
    if p0_to_p1 == 0 or p1_to_p2 == 0:
        raise ValueError("A filtered event contains a zero copy-number change")
    if p0_to_p1 * p1_to_p2 > 0:
        return "Further same-direction change"
    return "Opposite-direction change"


def build_segments(parent: str, child: str, p1_events: pd.DataFrame,
                   p2_events: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for p1 in p1_events.itertuples(index=False):
        candidates = p2_events[
            (p2_events["chr"] == p1.chr)
            & (p2_events["end_pos"] > p1.start_pos)
            & (p2_events["start_pos"] < p1.end_pos)
        ]
        overlaps: list[dict[str, object]] = []
        for p2 in candidates.itertuples(index=False):
            start = max(int(p1.start_pos), int(p2.start_pos))
            end = min(int(p1.end_pos), int(p2.end_pos))
            if end <= start:
                continue
            overlaps.append({
                "start": start,
                "end": end,
                "category": change_category(float(p1.cn_a), float(p1.cn_b),
                                            float(p2.cn_a), float(p2.cn_b)),
                "P2_event_id": p2.event_id,
                "P2_pair_P1_CN": float(p2.cn_a),
                "P2_CN": float(p2.cn_b),
                "P1_CN_consistent": bool(float(p1.cn_b) == float(p2.cn_a)),
            })
        overlaps.sort(key=lambda row: (int(row["start"]), int(row["end"])))
        for left, right in zip(overlaps, overlaps[1:]):
            if int(right["start"]) < int(left["end"]):
                raise ValueError(f"Overlapping child events would double-count {parent}-{child}")

        cursor = int(p1.start_pos)
        for overlap in overlaps:
            start = int(overlap["start"])
            end = int(overlap["end"])
            if start > cursor:
                records.append(segment_record(parent, child, p1, cursor, start,
                                              "Exact inheritance"))
            records.append(segment_record(
                parent, child, p1, start, end, str(overlap["category"]),
                p2_event_id=str(overlap["P2_event_id"]),
                p2_pair_p1_cn=float(overlap["P2_pair_P1_CN"]),
                p2_cn=float(overlap["P2_CN"]),
                p1_cn_consistent=bool(overlap["P1_CN_consistent"]),
            ))
            cursor = end
        if cursor < int(p1.end_pos):
            records.append(segment_record(parent, child, p1, cursor, int(p1.end_pos),
                                          "Exact inheritance"))

    result = pd.DataFrame(records)
    if result.empty:
        raise ValueError(f"No segments generated for {parent}-{child}")
    expected = int((p1_events["end_pos"] - p1_events["start_pos"]).sum())
    observed = int(result["length_bp"].sum())
    if observed != expected:
        raise RuntimeError(f"Segment lengths do not reconcile for {parent}-{child}: {observed} != {expected}")
    return result


def segment_record(parent: str, child: str, p1: object, start: int, end: int,
                   category: str, p2_event_id: str = "", p2_pair_p1_cn: float = np.nan,
                   p2_cn: float = np.nan, p1_cn_consistent: object = "") -> dict[str, object]:
    return {
        "parent": parent,
        "child": child,
        "P1_event_id": p1.event_id,
        "chr": p1.chr,
        "P1_event_start": int(p1.start_pos),
        "P1_event_end": int(p1.end_pos),
        "segment_start": start,
        "segment_end": end,
        "length_bp": end - start,
        "U1_CN": float(p1.cn_a),
        "P1_CN": float(p1.cn_b),
        "P2_event_id": p2_event_id,
        "P2_pair_P1_CN": p2_pair_p1_cn,
        "P2_CN": p2_cn,
        "P1_CN_consistent": p1_cn_consistent,
        "category": category,
    }


def summarize(details: pd.DataFrame, p1_event_counts: dict[str, int],
              p2_event_counts: dict[tuple[str, str], int]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for parent, children in BRANCHES.items():
        for child in children:
            group = details[(details["parent"] == parent) & (details["child"] == child)]
            denominator = int(group["length_bp"].sum())
            for category in CATEGORY_ORDER:
                category_bp = int(group.loc[group["category"] == category, "length_bp"].sum())
                rows.append({
                    "parent": parent,
                    "child": child,
                    "category": category,
                    "filtered_P1_CNV_events": p1_event_counts[parent],
                    "filtered_P1_P2_CNV_events": p2_event_counts[(parent, child)],
                    "P1_CNV_denominator_bp": denominator,
                    "category_bp": category_bp,
                    "bp_fraction": category_bp / denominator,
                    "percent": 100 * category_bp / denominator,
                })
    summary = pd.DataFrame(rows)
    totals = summary.groupby(["parent", "child"], sort=False)["bp_fraction"].sum()
    if not np.allclose(totals.to_numpy(), 1.0, atol=1e-12):
        raise RuntimeError("The three fate proportions do not sum to one")
    sibling_denominators = summary.groupby(["parent", "child"], sort=False)["P1_CNV_denominator_bp"].first()
    for parent in BRANCHES:
        values = sibling_denominators.loc[parent].to_numpy()
        if len(set(values.tolist())) != 1:
            raise RuntimeError(f"Sibling denominators differ for {parent}")
    return summary


def build_qc(details: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for parent, children in BRANCHES.items():
        for child in children:
            group = details[(details["parent"] == parent) & (details["child"] == child)]
            changed = group[group["category"] != "Exact inheritance"]
            mismatch = changed[changed["P1_CN_consistent"] == False]  # noqa: E712
            rows.append({
                "parent": parent,
                "child": child,
                "changed_overlap_bp": int(changed["length_bp"].sum()),
                "P1_CN_mismatch_bp": int(mismatch["length_bp"].sum()),
                "P1_CN_mismatch_fraction_of_denominator": (
                    float(mismatch["length_bp"].sum()) / float(group["length_bp"].sum())
                ),
            })
    return pd.DataFrame(rows)


def plot_stacked(summary: pd.DataFrame, output_dir: Path) -> None:
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 6.5,
        "axes.labelsize": 7.5,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(165 / 25.4, 85 / 25.4))
    positions = np.asarray([0, 1, 2.8, 3.8, 5.6, 6.6])
    pairs = [(parent, child) for parent in BRANCHES for child in BRANCHES[parent]]
    bottoms = np.zeros(len(pairs))

    for category in CATEGORY_ORDER:
        values = np.asarray([
            float(summary.loc[
                (summary["parent"] == parent)
                & (summary["child"] == child)
                & (summary["category"] == category),
                "bp_fraction",
            ].iloc[0])
            for parent, child in pairs
        ])
        bars = ax.bar(positions, values, bottom=bottoms, width=0.68,
                      color=COLORS[category], edgecolor="white", linewidth=0.5)
        for bar, value, bottom in zip(bars, values, bottoms):
            if value >= 0.035:
                ax.text(bar.get_x() + bar.get_width() / 2, bottom + value / 2,
                        f"{value * 100:.1f}%", ha="center", va="center", fontsize=6.0,
                        color="black" if category == "Further same-direction change" else "white")
        bottoms += values

    ax.set_xticks(positions, [child for _, child in pairs])
    for center, parent in zip((0.5, 3.3, 6.1), BRANCHES):
        ax.text(center, -0.145, parent, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.5, fontweight="bold", clip_on=False)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1, decimals=0))
    ax.set_ylabel("Proportion of high-confidence P1 CNV length", labelpad=11)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(width=0.5, length=2.2)
    handles = [Patch(facecolor=COLORS[name], label=name) for name in CATEGORY_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.015), handlelength=1.15,
               columnspacing=1.4, handletextpad=0.45)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.97, bottom=0.275)

    for suffix, kwargs in (("pdf", {}), ("png", {"dpi": 600})):
        path = output_dir / f"{OUTPUT_STEM}_stacked_bar.{suffix}"
        fig.savefig(path, format=suffix, facecolor="white", **kwargs)
        print(f"Saved {path}")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path,
                        help="Directory containing CNVresult_filter.csv")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    cnv_path = args.input / "CNVresult_filter.csv"
    cnv = pd.read_csv(cnv_path)
    required = ["pair_sampleA_name", "pair_sampleB_name", "chr", "start_pos",
                "end_pos", "length", "copy_number_A", "copy_number_B"]
    missing = [column for column in required if column not in cnv.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    cnv = cnv[cnv["chr"].isin([f"chr{i}" for i in range(1, 23)])].copy()

    detail_frames = []
    p1_event_counts: dict[str, int] = {}
    p2_event_counts: dict[tuple[str, str], int] = {}
    for parent, children in BRANCHES.items():
        p1_events = oriented_pair_events(cnv, "U1", parent)
        p1_event_counts[parent] = len(p1_events)
        for child in children:
            p2_events = oriented_pair_events(cnv, parent, child)
            p2_event_counts[(parent, child)] = len(p2_events)
            detail_frames.append(build_segments(parent, child, p1_events, p2_events))

    details = pd.concat(detail_frames, ignore_index=True)
    summary = summarize(details, p1_event_counts, p2_event_counts)
    qc = build_qc(details)
    summary_wide = summary.pivot(
        index=["parent", "child", "filtered_P1_CNV_events", "filtered_P1_P2_CNV_events",
               "P1_CNV_denominator_bp"],
        columns="category", values="bp_fraction",
    ).reset_index()
    summary_wide.columns.name = None

    write_csv(details, args.output / f"{OUTPUT_STEM}_interval_details.csv")
    write_csv(summary, args.output / f"{OUTPUT_STEM}_long.csv")
    write_csv(summary_wide, args.output / f"{OUTPUT_STEM}_summary.csv")
    write_csv(qc, args.output / f"{OUTPUT_STEM}_QC.csv")
    plot_stacked(summary, args.output)

    method = (
        "Analysis universe: filtered high-confidence U1-to-P1 CNV events in CNVresult_filter.csv.\n"
        "Denominator: the total union length of filtered U1-to-P1 CNV events for each P1 branch.\n"
        "Coordinates are treated as half-open [start_pos, end_pos), so length=end_pos-start_pos, matching the source length column.\n"
        "Exact inheritance (operational definition): P1-CNV bases not overlapped by any filtered P1-to-P2 CNV event.\n"
        "Further same-direction change: overlap length where the U1-to-P1 and P1-to-P2 event directions have the same sign.\n"
        "Opposite-direction change: overlap length where the two event directions have opposite signs.\n"
        "Filtered P1-to-P2 events outside the filtered U1-to-P1 CNV universe are excluded from the denominator.\n"
        "P1 CN state disagreements between the two pairwise calls are retained in the requested three-category result and flagged in QC.\n"
        "The three length-weighted fractions sum to 100% for each P1-P2 pair.\n"
    )
    (args.output / f"{OUTPUT_STEM}_METHOD.txt").write_text(method, encoding="utf-8")

    display = summary_wide.copy()
    for category in CATEGORY_ORDER:
        display[category] *= 100
    print(display.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    print("\nQC")
    print(qc.to_string(index=False, float_format=lambda value: f"{value:.6f}"))


if __name__ == "__main__":
    main()

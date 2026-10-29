#!/usr/bin/env python3
"""Split filtered P1-P2 CNV length into ancestral-overlap and P2-specific regions."""

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
CATEGORY_ORDER = ("Inherited CNV", "De novo CNV")
COLORS = {"Inherited CNV": "#4E79A7", "De novo CNV": "#F28E2B"}
OUTPUT_STEM = "six_P2_de_novo_vs_inherited_CNV_length"


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def oriented_pair_events(cnv: pd.DataFrame, sample_a: str, sample_b: str) -> pd.DataFrame:
    source_a = SAMPLE_COLUMNS[sample_a]
    source_b = SAMPLE_COLUMNS[sample_b]
    forward = cnv[
        (cnv["pair_sampleA_name"] == source_a)
        & (cnv["pair_sampleB_name"] == source_b)
    ].copy()
    forward["cn_a"] = pd.to_numeric(forward["copy_number_A"], errors="coerce")
    forward["cn_b"] = pd.to_numeric(forward["copy_number_B"], errors="coerce")

    reverse = cnv[
        (cnv["pair_sampleA_name"] == source_b)
        & (cnv["pair_sampleB_name"] == source_a)
    ].copy()
    reverse["cn_a"] = pd.to_numeric(reverse["copy_number_B"], errors="coerce")
    reverse["cn_b"] = pd.to_numeric(reverse["copy_number_A"], errors="coerce")

    events = pd.concat([forward, reverse], ignore_index=True)
    if events.empty:
        raise ValueError(f"No filtered events found for {sample_a}-{sample_b}")
    for column in ("start_pos", "end_pos", "length"):
        events[column] = pd.to_numeric(events[column], errors="raise").astype(np.int64)
    if events[["cn_a", "cn_b"]].isna().any().any():
        raise ValueError(f"Missing copy-number states in {sample_a}-{sample_b}")
    calculated_length = events["end_pos"] - events["start_pos"]
    if (calculated_length <= 0).any():
        raise ValueError(f"Non-positive event length in {sample_a}-{sample_b}")
    if not (calculated_length == events["length"]).all():
        raise ValueError(f"Source length does not equal end_pos-start_pos in {sample_a}-{sample_b}")

    events = events.sort_values(["chr", "start_pos", "end_pos"]).reset_index(drop=True)
    events["event_id"] = [f"{sample_a}_{sample_b}_HC{i + 1}" for i in range(len(events))]
    for chrom, group in events.groupby("chr", sort=False):
        previous_end = None
        for event in group.itertuples(index=False):
            if previous_end is not None and int(event.start_pos) < previous_end:
                raise ValueError(f"Overlapping filtered events in {sample_a}-{sample_b} on {chrom}")
            previous_end = int(event.end_pos)
    return events


def segment_record(parent: str, child: str, p2_event: object, start: int, end: int,
                   category: str, p1_event_id: str = "", p1_u1_cn: float = np.nan,
                   p1_cn: float = np.nan) -> dict[str, object]:
    return {
        "parent": parent,
        "child": child,
        "P2_event_id": p2_event.event_id,
        "chr": p2_event.chr,
        "P2_event_start": int(p2_event.start_pos),
        "P2_event_end": int(p2_event.end_pos),
        "segment_start": start,
        "segment_end": end,
        "length_bp": end - start,
        "P1_pair_CN": float(p2_event.cn_a),
        "P2_CN": float(p2_event.cn_b),
        "overlapping_P1_event_id": p1_event_id,
        "U1_CN_in_P1_event": p1_u1_cn,
        "P1_CN_in_P1_event": p1_cn,
        "category": category,
    }


def build_segments(parent: str, child: str, p1_events: pd.DataFrame,
                   p2_events: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for p2 in p2_events.itertuples(index=False):
        candidates = p1_events[
            (p1_events["chr"] == p2.chr)
            & (p1_events["end_pos"] > p2.start_pos)
            & (p1_events["start_pos"] < p2.end_pos)
        ]
        overlaps = []
        for p1 in candidates.itertuples(index=False):
            start = max(int(p1.start_pos), int(p2.start_pos))
            end = min(int(p1.end_pos), int(p2.end_pos))
            if end > start:
                overlaps.append((start, end, p1))
        overlaps.sort(key=lambda item: (item[0], item[1]))
        for left, right in zip(overlaps, overlaps[1:]):
            if right[0] < left[1]:
                raise ValueError(f"Overlapping P1 events would double-count {parent}-{child}")

        cursor = int(p2.start_pos)
        for start, end, p1 in overlaps:
            if start > cursor:
                records.append(segment_record(parent, child, p2, cursor, start, "De novo CNV"))
            records.append(segment_record(
                parent, child, p2, start, end, "Inherited CNV",
                p1_event_id=p1.event_id,
                p1_u1_cn=float(p1.cn_a),
                p1_cn=float(p1.cn_b),
            ))
            cursor = end
        if cursor < int(p2.end_pos):
            records.append(segment_record(parent, child, p2, cursor, int(p2.end_pos),
                                          "De novo CNV"))

    detail = pd.DataFrame(records)
    expected = int((p2_events["end_pos"] - p2_events["start_pos"]).sum())
    observed = int(detail["length_bp"].sum())
    if observed != expected:
        raise RuntimeError(f"P2 segment lengths do not reconcile for {parent}-{child}: {observed} != {expected}")
    return detail


def summarize(details: pd.DataFrame, p1_counts: dict[str, int],
              p2_counts: dict[tuple[str, str], int]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for parent, children in BRANCHES.items():
        for child in children:
            group = details[(details["parent"] == parent) & (details["child"] == child)]
            denominator = int(group["length_bp"].sum())
            if denominator <= 0:
                raise ValueError(f"No P2 CNV length for {parent}-{child}")
            for category in CATEGORY_ORDER:
                category_bp = int(group.loc[group["category"] == category, "length_bp"].sum())
                rows.append({
                    "parent": parent,
                    "child": child,
                    "category": category,
                    "filtered_U1_P1_CNV_events": p1_counts[parent],
                    "filtered_P1_P2_CNV_events": p2_counts[(parent, child)],
                    "P2_CNV_denominator_bp": denominator,
                    "category_bp": category_bp,
                    "bp_fraction": category_bp / denominator,
                    "percent": 100 * category_bp / denominator,
                })
    summary = pd.DataFrame(rows)
    totals = summary.groupby(["parent", "child"], sort=False)["bp_fraction"].sum().to_numpy()
    if not np.allclose(totals, 1.0, atol=1e-12):
        raise RuntimeError("Inherited and de novo proportions do not sum to one")
    return summary


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
                        f"{value * 100:.1f}%", ha="center", va="center",
                        fontsize=6.0, color="white")
        bottoms += values

    ax.set_xticks(positions, [child for _, child in pairs])
    for center, parent in zip((0.5, 3.3, 6.1), BRANCHES):
        ax.text(center, -0.145, parent, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.5, fontweight="bold", clip_on=False)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1, decimals=0))
    ax.set_ylabel("Proportion of P2 CNV length", labelpad=11)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(width=0.5, length=2.2)
    handles = [Patch(facecolor=COLORS[category], label=category) for category in CATEGORY_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.015), handlelength=1.15,
               columnspacing=1.5, handletextpad=0.45)
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

    cnv = pd.read_csv(args.input / "CNVresult_filter.csv")
    required = ["pair_sampleA_name", "pair_sampleB_name", "chr", "start_pos",
                "end_pos", "length", "copy_number_A", "copy_number_B"]
    missing = [column for column in required if column not in cnv.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    cnv = cnv[cnv["chr"].isin([f"chr{i}" for i in range(1, 23)])].copy()

    detail_frames = []
    p1_counts: dict[str, int] = {}
    p2_counts: dict[tuple[str, str], int] = {}
    for parent, children in BRANCHES.items():
        p1_events = oriented_pair_events(cnv, "U1", parent)
        p1_counts[parent] = len(p1_events)
        for child in children:
            p2_events = oriented_pair_events(cnv, parent, child)
            p2_counts[(parent, child)] = len(p2_events)
            detail_frames.append(build_segments(parent, child, p1_events, p2_events))

    details = pd.concat(detail_frames, ignore_index=True)
    summary = summarize(details, p1_counts, p2_counts)
    summary_wide = summary.pivot(
        index=["parent", "child", "filtered_U1_P1_CNV_events",
               "filtered_P1_P2_CNV_events", "P2_CNV_denominator_bp"],
        columns="category", values="bp_fraction",
    ).reset_index()
    summary_wide.columns.name = None

    write_csv(details, args.output / f"{OUTPUT_STEM}_interval_details.csv")
    write_csv(summary, args.output / f"{OUTPUT_STEM}_long.csv")
    write_csv(summary_wide, args.output / f"{OUTPUT_STEM}_summary.csv")
    plot_stacked(summary, args.output)

    method = (
        "Operational analysis based on filtered pairwise events in CNVresult_filter.csv.\n"
        "Denominator: total filtered P1-to-P2 CNV length for each P2 sample.\n"
        "Inherited CNV: P1-to-P2 CNV bases overlapping a filtered U1-to-P1 CNV interval.\n"
        "De novo CNV: P1-to-P2 CNV bases outside all filtered U1-to-P1 CNV intervals.\n"
        "Coordinates are half-open [start_pos, end_pos), matching source length=end_pos-start_pos.\n"
        "Overlapping and non-overlapping segments partition every P1-to-P2 event, so the two fractions sum to 100%.\n"
        "Important limitation: because no direct U1-to-P2 filtered comparisons are present, inherited means ancestral-interval overlap, not confirmed preservation of the same absolute CN state.\n"
    )
    (args.output / f"{OUTPUT_STEM}_METHOD.txt").write_text(method, encoding="utf-8")

    display = summary_wide.copy()
    for category in CATEGORY_ORDER:
        display[category] *= 100
    print(display.to_string(index=False, float_format=lambda value: f"{value:.4f}"))


if __name__ == "__main__":
    main()

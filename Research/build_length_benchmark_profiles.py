#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path


CHROMS = [f"chr{i}" for i in range(1, 23)]
LENGTHS = [20_000, 50_000, 100_000, 200_000, 500_000, 1_000_000]
TARGET_MEAN_COVERAGE = 30.0


@dataclass(frozen=True)
class Segment:
    chrom: str
    start: int
    end: int
    cn: float
    state: str
    source: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Event:
    event_id: str
    pair_id: str
    chrom: str
    start: int
    end: int
    a_cn: int
    b_cn: int

    @property
    def length(self) -> int:
        return self.end - self.start


def read_fasta_length(path: Path) -> int:
    length = 0
    header_seen = False
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header_seen:
                    raise ValueError(f"Expected one FASTA record in {path}")
                header_seen = True
            else:
                length += len(line)
    if not header_seen:
        raise ValueError(f"No FASTA header in {path}")
    return length


def load_source(path: Path, state: str, lengths: dict[str, int]) -> list[Segment]:
    result: list[Segment] = []
    with path.open() as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip().split("\t")
            if len(fields) != 4:
                raise ValueError(f"{path}:{line_no}: expected 4 columns")
            chrom, start_s, end_s, cn_s = fields
            if chrom not in lengths:
                continue
            start, end, cn = int(start_s), int(end_s), float(cn_s)
            if not (0 <= start < end <= lengths[chrom]):
                raise ValueError(f"Out-of-range interval {path}:{line_no}")
            result.append(Segment(chrom, start, end, cn, state, path.name))
    return result


def build_profile(calls: list[Segment], lengths: dict[str, int]) -> list[Segment]:
    by_chrom: dict[str, list[Segment]] = {chrom: [] for chrom in CHROMS}
    for call in calls:
        by_chrom[call.chrom].append(call)
    profile: list[Segment] = []
    for chrom in CHROMS:
        cursor = 0
        for call in sorted(by_chrom[chrom], key=lambda x: (x.start, x.end)):
            if call.start < cursor:
                raise ValueError(f"Overlapping source calls on {chrom} near {call.start}")
            if cursor < call.start:
                profile.append(Segment(chrom, cursor, call.start, 2.0, "NEUTRAL_ASSUMED", "gap"))
            profile.append(call)
            cursor = call.end
        if cursor < lengths[chrom]:
            profile.append(Segment(chrom, cursor, lengths[chrom], 2.0, "NEUTRAL_ASSUMED", "gap"))
    return profile


def load_master_events(path: Path) -> dict[int, list[Event]]:
    result = {length: [] for length in LENGTHS}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            length = int(row["length_bp"])
            if length not in result:
                continue
            event = Event(
                row["event_id"], row["pair_id"], row["chrom"],
                int(row["start"]), int(row["end"]),
                int(row["Sample_A_CN"]), int(row["Sample_B_CN"]),
            )
            if event.length != length:
                raise ValueError(f"Length mismatch for {event.event_id}")
            result[length].append(event)
    for length, events in result.items():
        dup = sum(event.b_cn > event.a_cn for event in events)
        dele = sum(event.b_cn < event.a_cn for event in events)
        if len(events) != 20 or dup != 10 or dele != 10:
            raise ValueError(f"Expected 10 dup and 10 del for {length}, got {dup} dup/{dele} del")
    return result


def apply_events(profile: list[Segment], events: list[Event]) -> list[Segment]:
    by_chrom: dict[str, list[Event]] = {chrom: [] for chrom in CHROMS}
    for event in events:
        by_chrom[event.chrom].append(event)
    output: list[Segment] = []
    applied: set[str] = set()
    for seg in profile:
        relevant = [e for e in by_chrom[seg.chrom] if seg.start <= e.start and e.end <= seg.end]
        cursor = seg.start
        for event in sorted(relevant, key=lambda x: x.start):
            if abs(seg.cn - event.a_cn) > 1e-9:
                raise ValueError(f"Baseline CN mismatch for {event.event_id}: {seg.cn} != {event.a_cn}")
            if cursor < event.start:
                output.append(replace(seg, start=cursor, end=event.start))
            output.append(Segment(seg.chrom, event.start, event.end, float(event.b_cn), "BENCHMARK", event.event_id))
            cursor = event.end
            applied.add(event.event_id)
        if cursor < seg.end:
            output.append(replace(seg, start=cursor))
    expected = {event.event_id for event in events}
    if applied != expected:
        raise ValueError(f"Unapplied events: {sorted(expected - applied)}")
    return output


def weighted_ploidy(profile: list[Segment]) -> float:
    total_bp = sum(seg.length for seg in profile)
    return sum(seg.length * seg.cn for seg in profile) / total_bp


def write_profile(path: Path, profile: list[Segment], ploidy: float) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "length_bp", "absolute_CN", "state", "source", "target_coverage"])
        for seg in profile:
            coverage = TARGET_MEAN_COVERAGE * seg.cn / ploidy
            writer.writerow([seg.chrom, seg.start, seg.end, seg.length, f"{seg.cn:g}", seg.state, seg.source, f"{coverage:.8f}"])


def write_benchmark(path: Path, events: list[Event]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "event_id", "pair_id", "chrom", "start", "end", "length_bp",
            "length_group", "Sample_A_CN", "Sample_B_CN", "delta_B_minus_A",
            "direction", "log2_CN_ratio", "abs_log2_CN_ratio",
        ])
        for event in sorted(events, key=lambda x: (CHROMS.index(x.chrom), x.start)):
            ratio = math.log2(event.b_cn / event.a_cn)
            writer.writerow([
                event.event_id, event.pair_id, event.chrom, event.start, event.end,
                event.length, f"{event.length // 1000}kb", event.a_cn, event.b_cn,
                event.b_cn - event.a_cn, "dup" if event.b_cn > event.a_cn else "del",
                f"{ratio:.12f}", f"{abs(ratio):.12f}",
            ])


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit("Usage: prepare.py GAIN LOSS REF_DIR MASTER_BENCHMARK OUT_ROOT BASELINE_PROFILE")
    gain_path, loss_path, ref_dir, master_path, out_root, baseline_reference = map(Path, sys.argv[1:])
    lengths = {chrom: read_fasta_length(ref_dir / f"{chrom}.fa") for chrom in CHROMS}
    calls = load_source(gain_path, "GAIN_SOURCE", lengths) + load_source(loss_path, "LOSS_SOURCE", lengths)
    profile_a = build_profile(calls, lengths)
    events_by_length = load_master_events(master_path)
    ploidy_a = weighted_ploidy(profile_a)

    for length in LENGTHS:
        label = f"{length // 1000}kb"
        experiment = out_root / "experiments" / label
        inputs = experiment / "inputs"
        benchmark = experiment / "benchmark"
        inputs.mkdir(parents=True, exist_ok=True)
        benchmark.mkdir(parents=True, exist_ok=True)
        events = events_by_length[length]
        profile_b = apply_events(profile_a, events)
        ploidy_b = weighted_ploidy(profile_b)
        if abs(ploidy_a - ploidy_b) > 1e-12:
            raise AssertionError(f"Ploidy mismatch for {label}: {ploidy_a} vs {ploidy_b}")
        write_profile(inputs / "Sample_A.absolute_CN_segments.tsv", profile_a, ploidy_a)
        write_profile(inputs / "Sample_B.absolute_CN_segments.tsv", profile_b, ploidy_b)
        write_benchmark(benchmark / "Sample_B_vs_Sample_A.benchmark.tsv", events)
        with (benchmark / "design_summary.tsv").open("w") as handle:
            handle.write("key\tvalue\n")
            handle.write(f"experiment_length_bp\t{length}\n")
            handle.write(f"length_group\t{label}\n")
            handle.write("benchmark_events\t20\n")
            handle.write("dup_events\t10\n")
            handle.write("del_events\t10\n")
            handle.write("event_CN_transition\tCN2_to_CN3_and_CN3_to_CN2\n")
            handle.write(f"Sample_A_length_weighted_ploidy\t{ploidy_a:.12f}\n")
            handle.write(f"Sample_B_length_weighted_ploidy\t{ploidy_b:.12f}\n")
            handle.write(f"ploidy_difference\t{ploidy_b-ploidy_a:.12g}\n")
            handle.write("coordinate_system\t0-based_half-open\n")
            handle.write("event_source\toriginal_combined_length_benchmark_filtered_by_length\n")
            handle.write(f"master_benchmark\t{master_path}\n")
            handle.write("unlisted_source_regions\tassumed_CN2\n")
            handle.write("coverage_formula\t30_times_segment_CN_divided_by_length_weighted_ploidy\n")

    generated_baseline = out_root / "experiments" / "20kb" / "inputs" / "Sample_A.absolute_CN_segments.tsv"
    if baseline_reference.exists() and generated_baseline.read_bytes() != baseline_reference.read_bytes():
        raise AssertionError("Generated Sample_A baseline differs from the reference project")


if __name__ == "__main__":
    main()

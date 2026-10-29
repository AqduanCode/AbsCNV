#!/usr/bin/env python3
from __future__ import annotations

import csv
import random
import sys
from dataclasses import dataclass, replace
from pathlib import Path


CHROMS = [f"chr{i}" for i in range(1, 23)]
TARGET_MEAN_COVERAGE = 30.0
SEED = 20260803
MIN_EVENT_GAP = 250_000
EDGE_MARGIN = 20_000
EVENT_COUNT = 50
MIN_EVENT_LENGTH = 20_000
MAX_EVENT_LENGTH = 1_000_000
EVENT_LENGTH_STEP = 10_000


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


def read_fasta(path: Path) -> tuple[str, str]:
    header = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"Expected one FASTA record in {path}")
                header = line[1:].split()[0]
            else:
                chunks.append(line.upper())
    if header is None:
        raise ValueError(f"No FASTA header in {path}")
    return header, "".join(chunks)


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


def overlaps_selected(chrom: str, start: int, end: int, selected: list[Event]) -> bool:
    for event in selected:
        if event.chrom != chrom:
            continue
        if start < event.end + MIN_EVENT_GAP and end + MIN_EVENT_GAP > event.start:
            return True
    return False


def choose_interval(
    cn: int,
    length: int,
    profile: list[Segment],
    sequences: dict[str, str],
    selected: list[Event],
    rng: random.Random,
) -> tuple[str, int, int] | None:
    eligible = [
        seg for seg in profile
        if abs(seg.cn - cn) <= 1e-9 and seg.length >= length + 2 * EDGE_MARGIN
    ]
    rng.shuffle(eligible)
    for seg in eligible:
        left = ((seg.start + EDGE_MARGIN + 9_999) // 10_000) * 10_000
        right = ((seg.end - EDGE_MARGIN - length) // 10_000) * 10_000
        if right < left:
            continue
        slot_count = (right - left) // 10_000 + 1
        slots = list(range(slot_count)) if slot_count <= 200 else rng.sample(range(slot_count), 200)
        rng.shuffle(slots)
        for slot in slots:
            start = left + slot * 10_000
            end = start + length
            if overlaps_selected(seg.chrom, start, end, selected):
                continue
            sequence = sequences[seg.chrom][start:end]
            if len(sequence) == length and set(sequence) <= set("ACGT"):
                return seg.chrom, start, end
    return None


def choose_events(profile: list[Segment], sequences: dict[str, str]) -> list[Event]:
    rng = random.Random(SEED)
    events: list[Event] = []
    for event_index in range(1, EVENT_COUNT + 1):
        chosen = None
        for _ in range(5000):
            a_cn = rng.randint(1, 6)
            b_cn = rng.choice([cn for cn in range(1, 7) if cn != a_cn])
            length = rng.randrange(
                MIN_EVENT_LENGTH // EVENT_LENGTH_STEP,
                MAX_EVENT_LENGTH // EVENT_LENGTH_STEP + 1,
            ) * EVENT_LENGTH_STEP
            region = choose_interval(a_cn, length, profile, sequences, events, rng)
            if region is not None:
                chosen = (a_cn, b_cn, length, region)
                break
        if chosen is None:
            raise RuntimeError(f"Could not place random event {event_index} of {EVENT_COUNT}")
        a_cn, b_cn, length, region = chosen
        design_id = f"R{event_index:02d}_CN{a_cn}_to_CN{b_cn}_{length//1000}kb"
        events.append(Event(f"E{event_index:02d}", design_id, *region, a_cn, b_cn))
    return sorted(events, key=lambda x: (CHROMS.index(x.chrom), x.start))


def apply_events(profile: list[Segment], events: list[Event]) -> list[Segment]:
    by_chrom: dict[str, list[Event]] = {chrom: [] for chrom in CHROMS}
    for event in events:
        by_chrom[event.chrom].append(event)
    output: list[Segment] = []
    for seg in profile:
        cursor = seg.start
        relevant = [e for e in by_chrom[seg.chrom] if seg.start <= e.start and e.end <= seg.end]
        for event in sorted(relevant, key=lambda x: x.start):
            if cursor < event.start:
                output.append(replace(seg, start=cursor, end=event.start))
            output.append(Segment(seg.chrom, event.start, event.end, float(event.b_cn), "BENCHMARK", event.event_id))
            cursor = event.end
        if cursor < seg.end:
            output.append(replace(seg, start=cursor))
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


def write_legacy_bed(path: Path, profile: list[Segment]) -> None:
    with path.open("w") as handle:
        handle.write("# chr\tchr_start\tchr_stop\tnum_positions\tcopy_number\n")
        for seg in profile:
            # This is an audit-format conversion only. Stock CNV-Sim must not be used
            # for loss segments because negative values alter the control side.
            variation = seg.cn - 2.0
            handle.write(f"{seg.chrom}\t{seg.start}\t{seg.end}\t{seg.length}\t{variation:g}\n")


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit("Usage: build_paired_benchmark_truth.py GAIN LOSS REF_DIR OUT_DIR")
    gain_path, loss_path, ref_dir, out_dir = map(Path, sys.argv[1:])
    out_dir.mkdir(parents=True, exist_ok=True)

    sequences: dict[str, str] = {}
    lengths: dict[str, int] = {}
    for chrom in CHROMS:
        header, sequence = read_fasta(ref_dir / f"{chrom}.fa")
        if header != chrom:
            raise ValueError(f"Expected {chrom}, found {header}")
        sequences[chrom] = sequence
        lengths[chrom] = len(sequence)

    calls = load_source(gain_path, "GAIN_SOURCE", lengths) + load_source(loss_path, "LOSS_SOURCE", lengths)
    profile_a = build_profile(calls, lengths)
    events = choose_events(profile_a, sequences)
    profile_b = apply_events(profile_a, events)
    ploidy_a = weighted_ploidy(profile_a)
    ploidy_b = weighted_ploidy(profile_b)
    write_profile(out_dir / "Sample_A.absolute_CN_segments.tsv", profile_a, ploidy_a)
    write_profile(out_dir / "Sample_B.absolute_CN_segments.tsv", profile_b, ploidy_b)
    write_legacy_bed(out_dir / "Sample_A.cnvsim_legacy.bed", profile_a)
    write_legacy_bed(out_dir / "Sample_B.cnvsim_legacy.bed", profile_b)

    with (out_dir / "Sample_B_vs_Sample_A.benchmark.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["event_id", "pair_id", "chrom", "start", "end", "length_bp", "Sample_A_CN", "Sample_B_CN", "delta_B_minus_A"])
        for e in events:
            writer.writerow([e.event_id, e.pair_id, e.chrom, e.start, e.end, e.length, e.a_cn, e.b_cn, e.b_cn - e.a_cn])

    with (out_dir / "design_summary.tsv").open("w") as handle:
        handle.write("key\tvalue\n")
        handle.write(f"target_mean_coverage\t{TARGET_MEAN_COVERAGE}\n")
        handle.write(f"Sample_A_length_weighted_ploidy\t{ploidy_a:.12f}\n")
        handle.write(f"Sample_B_length_weighted_ploidy\t{ploidy_b:.12f}\n")
        handle.write(f"ploidy_difference\t{ploidy_b-ploidy_a:.12g}\n")
        handle.write(f"benchmark_events\t{len(events)}\n")
        handle.write("event_design\t50_independent_random_CN_changes\n")
        handle.write("event_CN_range\t1_to_6\n")
        handle.write("event_length_range_bp\t20000_to_1000000_in_10000bp_steps\n")
        handle.write(f"random_seed\t{SEED}\n")
        handle.write("coordinate_system\t0-based_half-open\n")
        handle.write("unlisted_source_regions\tassumed_CN2\n")
        handle.write("coverage_formula\t30_times_segment_CN_divided_by_length_weighted_ploidy\n")


if __name__ == "__main__":
    main()

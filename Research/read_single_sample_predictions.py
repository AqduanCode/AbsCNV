"""Unmerged CNV readers with the manuscript's original coordinate conventions."""
from __future__ import annotations
import csv,math
from pathlib import Path
from dataclasses import dataclass

@dataclass(frozen=True)
class Call:
    call_id: str
    chrom: str
    start: int
    end: int
    direction: str
    source: str

    @property
    def length(self) -> int:
        return self.end - self.start

def normalize_chrom(value: str) -> str:
    token = value.strip().strip('"')
    if token.lower().startswith("chr"):
        token = token[3:]
    token = token.upper()
    if token == "23":
        token = "X"
    elif token == "24":
        token = "Y"
    return f"chr{token}"

def checked(call: Call) -> Call:
    if call.direction not in {"gain", "loss"}:
        raise ValueError(f"Invalid direction: {call}")
    if call.start < 0 or call.end <= call.start:
        raise ValueError(f"Invalid interval: {call}")
    return call

def read_event_tsv(path: Path) -> list[Call]:
    calls: list[Call] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"id", "chrom", "start", "end", "direction"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing required columns in {path}: {reader.fieldnames}")
        for row_no, row in enumerate(reader, 2):
            calls.append(checked(Call(
                row["id"].strip(), normalize_chrom(row["chrom"]),
                int(float(row["start"])), int(float(row["end"])),
                row["direction"].strip().lower(), f"{path}:{row_no}",
            )))
    if not calls:
        raise ValueError(f"No CNV events in {path}")
    return calls

def cn_direction(cn: float) -> str:
    if cn > 2:
        return "gain"
    if cn < 2:
        return "loss"
    return "neutral"

def read_mymethod(path: Path) -> list[Call]:
    calls = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chr", "start", "end", "Sample_A"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing MyMethod columns: {reader.fieldnames}")
        for row_no, row in enumerate(reader, 2):
            cn = float(row["Sample_A"])
            direc = cn_direction(cn)
            if direc == "neutral":
                raise ValueError(f"Neutral MyMethod row at {path}:{row_no}")
            calls.append(checked(Call(
                f"MyMethod_{row_no - 1}", normalize_chrom(row["chr"]),
                int(float(row["start"])) - 1, int(float(row["end"])), direc,
                f"{path}:{row_no};CN={cn:g}",
            )))
    return calls

def read_cnvnator(path: Path) -> list[Call]:
    calls = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"type", "chrom", "start", "end", "normalized_RD"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing CNVnator columns: {reader.fieldnames}")
        if "copy_number" in reader.fieldnames:
            raise ValueError("CNVnator raw input must not contain derived copy_number")
        for row_no, row in enumerate(reader, 2):
            kind = row["type"].lower()
            if kind not in {"deletion", "duplication"}:
                raise ValueError(f"Unexpected CNVnator type at {path}:{row_no}: {kind}")
            calls.append(checked(Call(
                f"CNVnator_{row_no - 1}", normalize_chrom(row["chrom"]),
                int(row["start"]) - 1, int(row["end"]),
                "loss" if kind == "deletion" else "gain",
                f"{path}:{row_no};normalized_RD={row['normalized_RD']}",
            )))
    return calls

def read_freec(path: Path) -> list[Call]:
    calls = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "start", "end", "copy_number", "status"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing ControlFREEC columns: {reader.fieldnames}")
        for row_no, row in enumerate(reader, 2):
            cn = float(row["copy_number"])
            direc = cn_direction(cn)
            status = row["status"].lower()
            if direc == "neutral" or status != direc:
                raise ValueError(f"ControlFREEC CN/status mismatch at {path}:{row_no}")
            calls.append(checked(Call(
                f"ControlFREEC_{row_no - 1}", normalize_chrom(row["chrom"]),
                int(row["start"]), int(row["end"]), direc,
                f"{path}:{row_no};CN={cn:g}",
            )))
    return calls

def read_cony(path: Path) -> list[Call]:
    calls = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "start", "end", "CN"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing CONY columns: {reader.fieldnames}")
        for row_no, row in enumerate(reader, 2):
            cn = float(row["CN"])
            direc = cn_direction(cn)
            if direc == "neutral":
                raise ValueError(f"Neutral CONY row at {path}:{row_no}")
            calls.append(checked(Call(
                f"CONY_{row_no - 1}", normalize_chrom(row["chrom"]),
                int(row["start"]) - 1, int(row["end"]), direc,
                f"{path}:{row_no};CN={cn:g}",
            )))
    return calls

def read_rdxplorer(path: Path) -> list[Call]:
    calls = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"state", "copyEst", "chrom", "posStart", "posEnd"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing RDXplorer columns: {reader.fieldnames}")
        for row_no, row in enumerate(reader, 2):
            state = int(float(row["state"]))
            copy_est = float(row["copyEst"])
            state_direction = "loss" if state == 1 else "gain" if state == 3 else "neutral"
            cn_direc = cn_direction(copy_est)
            if state_direction == "neutral" or cn_direc != state_direction:
                raise ValueError(f"RDXplorer state/copyEst mismatch at {path}:{row_no}")
            calls.append(checked(Call(
                f"RDXplorer_{row_no - 1}", normalize_chrom(row["chrom"]),
                int(float(row["posStart"])) - 1, int(float(row["posEnd"])), state_direction,
                f"{path}:{row_no};state={state};copyEst={copy_est:g}",
            )))
    return calls

#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import gzip
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ART = Path(os.environ.get("ABSCNV_ART", "/home/data3/home/jjy/CNVsim/CNV-Sim-py3/cnvsim/ART/art_illumina"))
READ_LENGTH = 100
FRAGMENT_MEAN = 200
FRAGMENT_SD = 10
MASTER_SEED = 20260923


def art_seed(sample: str, chrom: str, start: int, end: int, index: int) -> int:
    material = f"{MASTER_SEED}|{sample}|{chrom}|{start}|{end}|{index}".encode()
    return int.from_bytes(hashlib.sha256(material).digest()[:4], "big") % 2_147_483_646 + 1


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
                    raise ValueError(f"Expected one record in {path}")
                header = line[1:].split()[0]
            else:
                chunks.append(line.upper())
    if header is None:
        raise ValueError(f"No FASTA record in {path}")
    return header, "".join(chunks)


def append_gzip_member(source: Path, destination) -> int:
    lines = 0
    with source.open("rb") as src, gzip.GzipFile(fileobj=destination, mode="wb", compresslevel=6) as gz:
        while True:
            block = src.read(1024 * 1024)
            if not block:
                break
            lines += block.count(b"\n")
            gz.write(block)
    if lines % 4:
        raise ValueError(f"Malformed FASTQ {source}: {lines} lines")
    return lines // 4


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit("Usage: simulate_length_benchmark_reads.py SAMPLE CHROM SEGMENTS REF_FASTA OUT_DIR LOG")
    sample, chrom, segments_path, ref_path, out_dir, log_path = sys.argv[1:]
    segments_path, ref_path, out_dir, log_path = map(Path, [segments_path, ref_path, out_dir, log_path])
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    header, sequence = read_fasta(ref_path)
    if header != chrom:
        raise ValueError(f"Expected {chrom}, found {header}")

    rows = []
    with segments_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["chrom"] == chrom:
                rows.append(row)

    r1_path = out_dir / f"{sample}.{chrom}.R1.fastq.gz"
    r2_path = out_dir / f"{sample}.{chrom}.R2.fastq.gz"
    total_pairs = 0
    with r1_path.open("wb") as r1_out, r2_path.open("wb") as r2_out, log_path.open("w") as log:
        log.write("segment\tstart\tend\tCN\ttarget_coverage\tread_pairs\n")
        with tempfile.TemporaryDirectory(prefix=f"{sample}_{chrom}_", dir=str(out_dir)) as tmp_name:
            tmp = Path(tmp_name)
            for index, row in enumerate(rows, 1):
                start, end = int(row["start"]), int(row["end"])
                coverage = float(row["target_coverage"])
                segment = sequence[start:end]
                if not segment:
                    raise ValueError(f"Empty segment {chrom}:{start}-{end}")
                if coverage <= 0:
                    log.write(f"{index}\t{start}\t{end}\t{row['absolute_CN']}\t{coverage:.8f}\t0\n")
                    continue
                fasta = tmp / "segment.fa"
                prefix = tmp / "art_"
                with fasta.open("w") as handle:
                    handle.write(f">{sample}_{chrom}_{start}_{end}_seg{index}\n")
                    for offset in range(0, len(segment), 80):
                        handle.write(segment[offset:offset+80] + "\n")
                command = [
                    str(ART), "-i", str(fasta), "-p", "-m", str(FRAGMENT_MEAN),
                    "-s", str(FRAGMENT_SD), "-l", str(READ_LENGTH), "-f",
                    f"{coverage:.8f}", "-rs", str(art_seed(sample, chrom, start, end, index)),
                    "-na", "-o", str(prefix),
                ]
                completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                if completed.returncode != 0:
                    raise RuntimeError(f"ART failed for {chrom}:{start}-{end}: {completed.stderr}")
                pairs1 = append_gzip_member(Path(str(prefix) + "1.fq"), r1_out)
                pairs2 = append_gzip_member(Path(str(prefix) + "2.fq"), r2_out)
                if pairs1 != pairs2:
                    raise ValueError(f"Mate count mismatch at {chrom}:{start}-{end}")
                total_pairs += pairs1
                log.write(f"{index}\t{start}\t{end}\t{row['absolute_CN']}\t{coverage:.8f}\t{pairs1}\n")
                for artifact in tmp.glob("art_*"):
                    artifact.unlink()
                fasta.unlink()
    print(f"{sample}\t{chrom}\t{len(rows)}\t{total_pairs}\t{r1_path}\t{r2_path}")


if __name__ == "__main__":
    main()

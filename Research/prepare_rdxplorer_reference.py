#!/usr/bin/env python3
import re
import argparse
import shutil
from pathlib import Path

import numpy as np

RDX = OFFICIAL = PATCHED = REFERENCE = None
WINDOW = 100


def copy_code():
    if PATCHED.exists():
        raise FileExistsError("Use a new RDXplorer output directory: " + str(PATCHED))
    PATCHED.mkdir(parents=True)
    for path in OFFICIAL.iterdir():
        if path.name == 'rlib':
            target = PATCHED / 'rlib'
            target.mkdir()
            for item in path.iterdir():
                if item.is_file():
                    shutil.copy2(item, target / item.name)
        elif path.is_file():
            shutil.copy2(path, PATCHED / path.name)


def replace_once(path, old, new):
    text = path.read_text()
    if old not in text:
        raise RuntimeError(f'Expected patch target not found in {path}: {old!r}')
    path.write_text(text.replace(old, new, 1))


def patch_r_code():
    ewt = PATCHED / 'rlib' / 'ewtMain.R'
    replace_once(ewt, 'lower <- 1.5; upper <- 2.5; sig <- 0.000001',
                 'lower <- baseCopy - 0.5; upper <- baseCopy + 0.5; sig <- 0.000001')
    replace_once(ewt, 'copyZstat[gapProbes,1] <- 2',
                 'copyZstat[gapProbes,1] <- baseCopy')
    replace_once(ewt,
                 'summary <- getSummary(apphome, nzd,nzd[which(tristate==0)],copyZstat[,1],ch_digit, baseCopy=2)',
                 'summary <- getSummary(apphome, nzd,nzd[which(tristate==0)],copyZstat[,1],ch_digit, baseCopy=baseCopy)')
    replace_once(ewt, 'plotit(wrkgdir, nzd, summary, chr, filter, dirsep)',
                 '# plot disabled in the hg38 batch port; tabular calls are retained')


def fasta_records(path):
    name = None
    chunks = []
    with path.open('rb') as handle:
        for line in handle:
            if line.startswith(b'>'):
                if name is not None:
                    yield name, b''.join(chunks).upper()
                name = line[1:].split(None, 1)[0].decode()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, b''.join(chunks).upper()


def write_gc_and_gaps():
    gc_dir = PATCHED / 'rlib' / f'gc{WINDOW}_hg38'
    gap_dir = PATCHED / 'rlib' / 'gap_hg38'
    gc_dir.mkdir(parents=True, exist_ok=True)
    gap_dir.mkdir(parents=True, exist_ok=True)
    expected = {f'chr{i}' for i in range(1, 23)}
    observed = set()
    size_rows = []
    for chrom, seq in fasta_records(REFERENCE):
        if chrom not in expected:
            continue
        observed.add(chrom)
        n_full = len(seq) // WINDOW
        matrix = np.frombuffer(seq[:n_full * WINDOW], dtype=np.uint8).reshape(n_full, WINDOW)
        gc = ((matrix == ord('G')) | (matrix == ord('C'))).sum(axis=1)
        np.savetxt(gc_dir / f'{chrom}gc{WINDOW}.txt', gc, fmt='%d')
        with (gap_dir / f'{chrom}_gap.txt').open('w') as out:
            idx = 0
            full_end = n_full * WINDOW
            for match in re.finditer(b'N+', seq):
                start, end = match.span()
                if start >= full_end:
                    continue
                end = min(end, full_end)
                idx += 1
                out.write(f'{idx}\t{chrom}\t{start}\t{end}\t{idx}\tN\t{end-start}\tN_run\tno\n')
        size_rows.append((chrom, len(seq), n_full))
    if observed != expected:
        raise RuntimeError(f'Reference chromosomes mismatch: missing={sorted(expected-observed)}')
    with (RDX / 'hg38_precompute.tsv').open('w') as out:
        out.write('chrom\tlength_bp\tfull_100bp_windows\n')
        for chrom, length, n_full in size_rows:
            out.write(f'{chrom}\t{length}\t{n_full}\n')


def write_patch_notes():
    (RDX / 'PATCHES.md').write_text(
        '# RDXplorer 3.2 GRCh38 compatibility port\n\n'
        '- The upstream RDXplorer 3.2 source is supplied explicitly with --official-source; retain its upstream provenance.\n'
        '- `patched/` contains only source and R code; hg18/hg19 precomputed arrays are not copied.\n'
        '- GC counts and N-gap intervals are regenerated from the exact GRCh38 FASTA used for alignment, at the official 100 bp window size.\n'
        '- Read counts are generated with modern `samtools depth` at mapping quality >=20, then aggregated to 100 bp windows. This uses the documented q20-only mode (`sumWithZero=False`) and avoids the removed `samtools pileup` command.\n'
        '- The sample baseline is set to baseCopy 2 by call_rdxplorer_chromosome.R. Hard-coded CN2 thresholds and gap state in `ewtMain.R` are replaced with `baseCopy`.\n'
        '- Plot generation is disabled; all tabular `.sum`, `.sum.filt`, `.ewt`, and normalized-depth outputs are retained.\n'
    )


def main():
    global RDX, OFFICIAL, PATCHED, REFERENCE
    parser = argparse.ArgumentParser()
    parser.add_argument("--official-source", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    RDX = args.output_dir; OFFICIAL = args.official_source
    REFERENCE = args.reference; PATCHED = RDX / "patched"
    RDX.mkdir(parents=True, exist_ok=True)
    if not OFFICIAL.is_dir() or not REFERENCE.is_file():
        raise FileNotFoundError('Official RDXplorer source or project reference is missing')
    copy_code()
    patch_r_code()
    write_gc_and_gaps()
    write_patch_notes()
    print(PATCHED)


if __name__ == '__main__':
    main()

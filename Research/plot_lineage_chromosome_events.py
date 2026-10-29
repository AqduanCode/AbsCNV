"""Plot positional CNV inheritance summaries for the three U1 branches.

Input: filtered, pairwise CNV calls from the AbsCNV self-test analysis.
Each top-matrix entry is the summed length of calls on that chromosome (Mb).
Lower tracks show the genomic positions of the original calls, not inherited
projections. A chromosome is enlarged when both the founder comparison and
at least one descendant comparison contain a call on it.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.font_manager import fontManager
from matplotlib.patches import ConnectionPatch, Rectangle
import numpy as np
import pandas as pd



BRANCHES = {
    "U1.1": ("G6A7first", ("G6A7E1", "G6A7E11")),
    "U1.2": ("G6B11first", ("G6B11C1", "G6B11C6")),
    "U1.3": ("G6C7first", ("G6C7E1", "G6C7E9")),
}
DISPLAY = {
    "G6": "U1",
    "G6A7first": "U1.1",
    "G6A7E1": "U1.1.1",
    "G6A7E11": "U1.1.2",
    "G6B11first": "U1.2",
    "G6B11C1": "U1.2.1",
    "G6B11C6": "U1.2.2",
    "G6C7first": "U1.3",
    "G6C7E1": "U1.3.1",
    "G6C7E9": "U1.3.2",
}
COLORS = ("#d57c3d", "#347fa2", "#45a58f")
DARK = "#40545c"
LIGHT = "#dfe5e7"
CHROMS = tuple(f"chr{i}" for i in range(1, 23))


def mm(value: float) -> float:
    return value / 25.4


def blend(color: str, amount: float) -> tuple[float, float, float]:
    rgb = np.array(to_rgb(color))
    return tuple((1 - amount) + amount * rgb)


def fmt_mb(value: float) -> str:
    if value >= 10:
        return f"{value:.1f}"
    if value >= 0.1:
        return f"{value:.1f}"
    return f"{value:.2f}"


def subset(data: pd.DataFrame, sample_a: str, sample_b: str) -> pd.DataFrame:
    out = data.loc[
        (data["pair_sampleA_name"] == sample_a)
        & (data["pair_sampleB_name"] == sample_b)
    ].copy()
    if out.empty:
        raise ValueError(f"No CNV calls for {sample_a} -> {sample_b}")
    return out


def add_axes_mm(fig, x: float, y_bottom: float, w: float, h: float, H: float):
    return fig.add_axes((x / 175, y_bottom / H, w / 175, h / H))


def draw_matrix(fig, tracks, labels, zoom, height_mm):
    # The matrix occupies exactly 143 mm, giving each chromosome a 6.5 mm cell.
    left, width, top, row_h = 30.0, 143.0, height_mm - 13.5, 8.2
    cell_w = width / 22
    bottom = top - 3 * row_h
    values = np.zeros((3, 22))
    for i, data in enumerate(tracks):
        for j, chrom in enumerate(CHROMS):
            values[i, j] = data.loc[data["chr"] == chrom, "length"].sum() / 1e6
    vmax = max(values.max(), 1)
    for j, chrom in enumerate(CHROMS):
        x = left + (j + 0.5) * cell_w
        fig.text(x / 175, (top + 2.2) / height_mm, chrom,
                 ha="center", va="bottom", fontsize=6.1, color="#202b2f")
        for i in range(3):
            y = top - (i + 1) * row_h + 0.55
            v = values[i, j]
            # Restore the original length-dependent intensity. Empty cells
            # remain near white; only positive CNV lengths use the gradient.
            shade = 0 if v <= 0 else 0.3 + 0.62 * np.log1p(v) / np.log1p(vmax)
            face = (0.992, 0.995, 0.996) if v <= 0 else blend(COLORS[i], shade)
            rect = Rectangle(((left + j * cell_w + 0.25) / 175, y / height_mm),
                             (cell_w - 0.5) / 175, (row_h - 1.1) / height_mm,
                             transform=fig.transFigure,
                             facecolor=face, edgecolor=LIGHT, linewidth=0.45)
            fig.add_artist(rect)
            if v > 0:
                fig.text(x / 175, (y + (row_h - 1.1) / 2) / height_mm,
                         fmt_mb(v), ha="center", va="center", fontsize=6.0,
                         fontweight="normal", color="black")
    for i, label in enumerate(labels):
        fig.text((left - 2) / 175, (top - (i + 0.5) * row_h) / height_mm,
                 label, ha="right", va="center", fontsize=7)
    for k, chrom in enumerate(zoom, start=1):
        j = CHROMS.index(chrom)
        rect = Rectangle(((left + j * cell_w + 0.08) / 175,
                          (bottom + 0.25) / height_mm),
                         (cell_w - 0.16) / 175, (3 * row_h - 0.5) / height_mm,
                         transform=fig.transFigure, fill=False, edgecolor=DARK,
                         linewidth=0.9)
        fig.add_artist(rect)
        if len(zoom) > 2:
            # Numbered callouts remain legible when four or seven chromosomes
            # are expanded; individual curved arrows would cross repeatedly.
            fig.text((left + (j + 0.5) * cell_w) / 175,
                     (bottom - 2.45) / height_mm, str(k), ha="center", va="center",
                     fontsize=6.2, color="white", fontweight="normal",
                     bbox={"boxstyle": "circle,pad=0.15", "fc": DARK, "ec": "none"})
    return left, cell_w, bottom, vmax


def draw_gradient_key(fig, labels, vmax, height_mm, matrix_bottom):
    # Position the compact key beneath the upper chromosome matrix, to the
    # right of both curved magnification arrows in the two-panel figure.
    fig.text(156.0 / 175, (matrix_bottom - 4.8) / height_mm,
             "CNV length (Mb)", ha="center", va="center", fontsize=5.0,
             color="#344147")
    x, width = 145.0, 23.0
    for i, (color, label) in enumerate(zip(COLORS, labels)):
        y = matrix_bottom - 9.0 - i * 4.0
        fig.text(142.0 / 175, (y + 1.15) / height_mm,
                 label, ha="right", va="center", fontsize=5.0)
        # Sample uniformly in log1p space, matching the matrix's color
        # mapping. Uniform sampling in Mb hid the pale half of the gradient.
        shade = 0.3 + 0.62 * np.linspace(0, 1, 256)
        gradient = np.array([blend(color, s) for s in shade])[None, :, :]
        ax = add_axes_mm(fig, x, y, width, 2.3, height_mm)
        ax.imshow(gradient, aspect="auto", extent=[0, 1, 0, 1],
                  origin="lower", interpolation="bilinear")
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#d9e0e3")
            spine.set_linewidth(0.5)
    fig.text(x / 175, (matrix_bottom - 20.9) / height_mm, ">0",
             ha="left", va="center", fontsize=5.0, color="#344147")
    fig.text((x + width / 2) / 175, (matrix_bottom - 20.9) / height_mm,
             f"{np.expm1(np.log1p(vmax) / 2):.1f}", ha="center", va="center",
             fontsize=5.0, color="#344147")
    fig.text((x + width) / 175, (matrix_bottom - 20.9) / height_mm,
             f"{vmax:.1f}", ha="right", va="center", fontsize=5.0,
             color="#344147")


def draw_track(fig, tracks, labels, chrom, max_pos, index, nzoom, height_mm,
               matrix_left, cell_w, matrix_bottom):
    row, col = divmod(index, 2)
    last_unpaired = nzoom % 2 == 1 and index == nzoom - 1
    x, w = ((68.5, 61.0) if last_unpaired else
            ((30.0, 61.0) if col == 0 else (107.0, 61.0)))
    y = height_mm - 62.0 - 43.0 * row - 28.0
    ax = add_axes_mm(fig, x, y, w, 28.0, height_mm)
    extent = max_pos[chrom] / 1e6
    ax.set_xlim(0, extent)
    ax.set_ylim(-0.55, 2.58)
    for i, data in enumerate(tracks):
        ypos = 2 - i
        ax.hlines(ypos, 0, extent, color=LIGHT, linewidth=2.2, zorder=1)
        for ev in data.loc[data["chr"] == chrom].itertuples():
            start, end = ev.start_pos / 1e6, ev.end_pos / 1e6
            ax.add_patch(Rectangle((start, ypos - 0.125), end - start, 0.25,
                                   facecolor=COLORS[i], edgecolor="none", zorder=3))
            if end - start < extent * 0.007:
                ax.vlines((start + end) / 2, ypos - 0.13, ypos + 0.13,
                          colors=COLORS[i], linewidth=1.25, zorder=4)
        if col == 0:
            ax.text(-0.018, ypos, labels[i], transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=6.1)
    ax.set_yticks([])
    for spine_name in ("left", "right", "top"):
        ax.spines[spine_name].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(axis="x", labelsize=6.1, length=2.6, pad=1.6)
    ax.grid(axis="x", color="#edf0f1", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"{chrom} (Mb)", fontsize=7, labelpad=3)
    if nzoom > 2:
        fig.text((x - 4.6) / 175, (y + 25.7) / height_mm, str(index + 1),
                 ha="center", va="center", fontsize=6.2, color="white",
                 fontweight="normal",
                 bbox={"boxstyle": "circle,pad=0.15", "fc": DARK, "ec": "none"})
    else:
        matrix_j = CHROMS.index(chrom)
        source_x = (matrix_left + (matrix_j + 0.5) * cell_w) / 175
        source_y = (matrix_bottom - 1.1) / height_mm
        dest_x = (x + w * 0.5) / 175
        dest_y = (y + 29.7) / height_mm
        arrow = ConnectionPatch((source_x, source_y), (dest_x, dest_y),
                                "figure fraction", "figure fraction",
                                axesA=fig, axesB=fig, arrowstyle="-|>",
                                mutation_scale=8, linewidth=0.85,
                                color="#8a9ba0", connectionstyle="arc3,rad=-0.10")
        fig.add_artist(arrow)


def plot_branch(data, max_pos, branch, outdir):
    ancestor, descendants = BRANCHES[branch]
    pair_names = (("G6", ancestor), (ancestor, descendants[0]),
                  (ancestor, descendants[1]))
    labels = [f"{DISPLAY[a]} → {DISPLAY[b]}" for a, b in pair_names]
    tracks = [subset(data, a, b) for a, b in pair_names]
    founder = set(tracks[0]["chr"])
    descendant = set(tracks[1]["chr"]) | set(tracks[2]["chr"])
    zoom = sorted(founder & descendant, key=lambda c: int(c[3:]))
    if not zoom:
        raise ValueError(f"No shared chromosomes for {branch}")
    nrows = math.ceil(len(zoom) / 2)
    height_mm = 64 + 43 * nrows + 7
    fig = plt.figure(figsize=(mm(175), mm(height_mm)), facecolor="white")
    left, cell_w, bottom, vmax = draw_matrix(fig, tracks, labels, zoom, height_mm)
    for i, chrom in enumerate(zoom):
        draw_track(fig, tracks, labels, chrom, max_pos, i, len(zoom),
                   height_mm, left, cell_w, bottom)
    draw_gradient_key(fig, labels, vmax, height_mm, bottom)
    stem = "CNV_inheritance_" + branch.replace(".", "")
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=600, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    print(f"{branch}: {len(zoom)} enlarged: {', '.join(zoom)}; "
          f"175 × {height_mm} mm; {png.name}, {pdf.name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cnv", type=Path, required=True)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).parent / "inheritance_figures_final")
    parser.add_argument("--font-dir", type=Path, default=Path(r"C:\Windows\Fonts"))
    args = parser.parse_args()
    for file in (args.font_dir / "arial.ttf", args.font_dir / "arialbd.ttf"):
        fontManager.addfont(str(file))
    plt.rcParams.update({"font.family": "Arial", "font.size": 7,
                         "pdf.fonttype": 42, "axes.unicode_minus": False})
    data = pd.read_csv(args.cnv)
    background = pd.read_csv(args.background, usecols=["chr", "end"])
    if (data["end_pos"] <= data["start_pos"]).any():
        raise ValueError("CNV calls with nonpositive intervals")
    if not np.allclose(data["length"], data["end_pos"] - data["start_pos"]):
        raise ValueError("CNV length does not match end-start")
    max_pos = background.groupby("chr")["end"].max().to_dict()
    for chrom in CHROMS:
        if chrom not in max_pos:
            raise ValueError(f"Missing chromosome extent: {chrom}")
        max_pos[chrom] = max(max_pos[chrom], data.loc[data.chr == chrom, "end_pos"].max())
    args.output.mkdir(parents=True, exist_ok=True)
    for branch in BRANCHES:
        plot_branch(data, max_pos, branch, args.output)


if __name__ == "__main__":
    main()

"""Plot Table 1 benchmark performance as paired F1 lollipop panels.

Run with Python 3 and matplotlib. The script reads the original evaluation
files, validates their metric consistency, and
writes a vector PDF and a 600-dpi PNG beside this script by default.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


# Data paths are set from --data-root or explicit metric-file arguments in main.
SINGLE = PAIRED = None
FONT_REGULAR = Path(r"C:\Windows\Fonts\arial.ttf")
FONT_BOLD = Path(r"C:\Windows\Fonts\arialbd.ttf")

TOOLS = {
    "Single sample": ["AbsCNV", "CNVnator", "Control-FREEC", "RDXplorer", "CONY"],
    "Paired sample": ["AbsCNV", "BICseq2", "CNVkit", "Control-FREEC", "ichorCNA"],
}

COLORS = {"Event": "#347A9C", "Base-pair": "#CF7939"}
FIG_WIDTH_MM = 175.0
FIG_HEIGHT_MM = 113.0
MIN_FONT_PT = 6.2


def read_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def load_metrics() -> dict[tuple[str, str, str], tuple[float, float, float]]:
    data = {}
    for design, source in (("Single sample", SINGLE), ("Paired sample", PAIRED)):
        for row in read_rows(source, ","):
            level = {"single80_many_to_many": "Event", "base_pair": "Base-pair"}.get(row["criterion"])
            if level is None:
                continue
            tool = row["tool"].replace("MyMethod", "AbsCNV").replace("AbsCNV_background_updated", "AbsCNV").replace("ControlFREEC", "Control-FREEC")
            data[(design, level, tool)] = tuple(float(row[m]) for m in ("precision", "recall", "f1"))

    expected_keys = {(sample, level, tool) for sample, tools in TOOLS.items() for level in ("Event", "Base-pair") for tool in tools}
    assert set(data) == expected_keys, f"Unexpected or missing source rows: {set(data) ^ expected_keys}"
    for key, values in data.items():
        p, r, f1 = values
        assert all(0 <= value <= 1 for value in values), key
        assert math.isclose(f1, 2 * p * r / (p + r), abs_tol=0.0005), key
    return data


def configure_fonts() -> None:
    for path in (FONT_REGULAR, FONT_BOLD):
        if not path.is_file():
            raise FileNotFoundError(f"Required Arial font missing: {path}")
        font_manager.fontManager.addfont(str(path))
    resolved = font_manager.findfont(font_manager.FontProperties(family="Arial"))
    assert Path(resolved).name.lower() == "arial.ttf", f"Arial not resolved: {resolved}"
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 7.2,
        "axes.titlesize": 8.2, "axes.labelsize": 7.2,
        "xtick.labelsize": 6.8, "ytick.labelsize": 7.0,
        "legend.fontsize": 7.1, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def plot(data: dict, output_dir: Path) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4), facecolor="white")
    for sample, left in (("Single sample", 0.14), ("Paired sample", 0.64)):
        ax = fig.add_axes([left, 0.125, 0.35, 0.735])
        ax.set_xlim(-0.06, 2.15)
        ax.set_ylim(9.85, -0.85)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1"])
        ax.set_yticks([2 * i + 0.5 for i in range(5)])
        ax.set_yticklabels(TOOLS[sample])
        ax.tick_params(axis="x", length=2.2, width=0.6, pad=2.5, labelsize=6.2)
        ax.tick_params(axis="y", length=0, pad=3.5, labelsize=6.6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_position(("data", 0))
        ax.spines["left"].set_linewidth(0.7)
        ax.spines["left"].set_color("#AAB5BB")
        ax.spines["bottom"].set_linewidth(0.7)
        ax.spines["bottom"].set_color("#555D63")
        ax.spines["bottom"].set_bounds(0, 1.24)
        ax.set_title(sample, loc="left", pad=7, fontsize=8.5, fontweight="bold")
        ax.text(1.45, -0.68, "Precision", ha="center", va="center", fontsize=6.5)
        ax.text(1.84, -0.68, "Recall", ha="center", va="center", fontsize=6.5)
        for label in ax.get_yticklabels():
            if label.get_text() == "AbsCNV":
                label.set_fontweight("bold")

        for i, tool in enumerate(TOOLS[sample]):
            for offset, level in ((0, "Event"), (1, "Base-pair")):
                y = 2 * i + offset
                precision, recall, f1 = data[(sample, level, tool)]
                color = COLORS[level]
                ax.plot([0, f1], [y, y], color=color, linewidth=0.72,
                        linestyle=(0, (2.2, 2.2)), zorder=1)
                ax.plot(f1, y, "o", markersize=5.8, color=color, zorder=3)
                ax.text(f1 + 0.05, y, f"{f1:.3f}", va="center", ha="left", fontsize=6.2,
                        color="#3E464A", bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.35})
                for center, value in ((1.45, precision), (1.84, recall)):
                    x0, x1 = (1.32, 1.58) if center == 1.45 else (1.71, 1.97)
                    ax.add_patch(Rectangle((x0, y - 0.34), x1 - x0, 0.68,
                                           facecolor="white", edgecolor="#AAB5BB", linewidth=0.7))
                    ax.text(center, y, f"{value:.3f}", ha="center", va="center", fontsize=6.2,
                            color="#3E464A")

    handles = [Line2D([0], [0], color=COLORS[level], linestyle=(0, (2.2, 2.2)),
                      linewidth=0.72, marker="o", markersize=5.8, label=f"{level} F1")
               for level in ("Event", "Base-pair")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.012),
               ncol=2, frameon=False, fontsize=6.8, columnspacing=2.2, handlelength=2.4)
    for text in fig.findobj(match=matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() >= MIN_FONT_PT, (text.get_text(), text.get_fontsize())
            assert text.get_fontfamily()[0] == "Arial", (text.get_text(), text.get_fontfamily())

    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "Figure2_benchmark_performance.pdf"
    png_path = output_dir / "Figure2_benchmark_performance.png"
    fig.savefig(pdf_path, format="pdf", bbox_inches=None)
    fig.savefig(png_path, format="png", dpi=600, bbox_inches=None)
    plt.close(fig)
    return pdf_path, png_path


def main() -> None:
    global SINGLE, PAIRED, FONT_REGULAR, FONT_BOLD
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, help="Extracted Zenodo datasets directory")
    parser.add_argument("--single-metrics", type=Path)
    parser.add_argument("--paired-metrics", type=Path)
    parser.add_argument("--font-dir", type=Path, default=Path(r"C:\Windows\Fonts"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = args.data_root
    if root is None:
        candidates = [repo / "Zenodo" / "datasets", repo.parent / "datasets", repo.parent / "Zenodo" / "datasets"]
        root = next((p for p in candidates if p.is_dir()), None)
    if root is None and not (args.single_metrics and args.paired_metrics):
        parser.error("Supply --data-root or both metric files")
    SINGLE = args.single_metrics or root / "single_sample/evaluation/single_sample_four_standard_summary.csv"
    PAIRED = args.paired_metrics or root / "paired_samples/evaluation/four_standard_f1_summary.csv"
    FONT_REGULAR = args.font_dir / "arial.ttf"
    FONT_BOLD = args.font_dir / "arialbd.ttf"
    configure_fonts()
    data = load_metrics()
    pdf_path, png_path = plot(data, args.output_dir)
    audit = args.output_dir / "Figure2_plotted_metrics.csv"
    with audit.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_design", "metric_level", "tool", "precision", "recall", "f1"])
        for key, values in data.items():
            writer.writerow([*key, *values])
    print(f"Validated {len(data)} source metric rows. No manuscript or old Table 1 is required.")
    print(f"Figure size: {FIG_WIDTH_MM:g} x {FIG_HEIGHT_MM:g} mm; all text >= {MIN_FONT_PT:g} pt; Arial.")
    print(pdf_path)
    print(png_path)


if __name__ == "__main__":
    main()

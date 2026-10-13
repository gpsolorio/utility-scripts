#!/usr/bin/env python3
"""Generate a read-length QC figure from a BAM file with summary cards."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator

NAVY = "#20344A"
BLUE = "#3977A8"
PALE = "#F0F4F7"
BORDER = "#DCE5ED"
GREY = "#64778E"
GRID = "#DCE5ED"
WHITE = "#FFFFFF"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a read-length QC figure from a BAM file."
    )
    parser.add_argument("--bam", required=True, type=Path, help="Input BAM file.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--filename", default="read_length_qc.png",
        help="Output filename (default: read_length_qc.png).",
    )
    args = parser.parse_args()
    if Path(args.filename).name != args.filename:
        parser.error("--filename must be a filename, not a path; use --outdir for the destination.")
    return args


def extract_read_lengths(bam_file):
    lengths = []
    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            length = read.query_length
            if length is not None and length > 0:
                lengths.append(length)
    if not lengths:
        raise ValueError("No mapped primary reads with valid query lengths were found.")
    return np.asarray(lengths, dtype=np.int64)


def calculate_statistics(lengths):
    sorted_lengths = np.sort(lengths)
    total_bases = int(sorted_lengths.sum())
    cumulative_bases = np.cumsum(sorted_lengths[::-1])
    n50_index = np.searchsorted(cumulative_bases, total_bases / 2)
    return {
        "n": int(len(sorted_lengths)),
        "mean": float(np.mean(sorted_lengths)),
        "median": float(np.median(sorted_lengths)),
        "n50": int(sorted_lengths[::-1][n50_index]),
        "total_bases": total_bases,
    }


def add_metric_card(ax, label, value, unit, precision=0):
    """Draw one summary card, formatting integer and decimal metrics safely."""
    ax.set_facecolor(PALE)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(BORDER)
        spine.set_linewidth(1.5)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    display_value = f"{value:,.{precision}f}" if precision else f"{int(value):,}"

    ax.text(
        0.12, 0.76, label.upper(), transform=ax.transAxes,
        fontsize=10, fontweight="bold", color=GREY, ha="left", va="center",
    )
    ax.text(
        0.12, 0.48, display_value, transform=ax.transAxes,
        fontsize=23, fontweight="bold", color=NAVY, ha="left", va="center",
    )
    ax.text(
        0.12, 0.20, unit, transform=ax.transAxes,
        fontsize=11, color=GREY, ha="left", va="center",
    )


def make_figure(lengths, stats, output_file):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.facecolor": WHITE,
        "figure.facecolor": WHITE,
        "savefig.facecolor": WHITE,
    })

    fig = plt.figure(figsize=(13, 7), facecolor=WHITE)
    outer = GridSpec(
        1, 2, figure=fig, width_ratios=[1.25, 1.0], wspace=0.12,
        left=0.07, right=0.97, top=0.83, bottom=0.14,
    )

    # Read-length histogram: bins start at 0 bp; no median line or legend.
    ax = fig.add_subplot(outer[0, 0])
    high = int(max(lengths))
    bin_width = 10
    bin_end = max(bin_width, ((high // bin_width) + 1) * bin_width)
    bins = np.arange(0, bin_end + bin_width, bin_width)

    ax.hist(lengths, bins=bins, color=BLUE, edgecolor=WHITE, linewidth=0.8)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Read length (bp)", fontsize=11, labelpad=9)
    ax.set_ylabel("Number of reads", fontsize=11, labelpad=9)
    ax.set_title("Read Length Distribution", fontsize=14, color=NAVY, pad=16)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    ax.tick_params(axis="both", colors="black")

    # Summary cards: three across the top, two centered below.
    summary_grid = outer[0, 1].subgridspec(2, 6, wspace=0.22, hspace=0.28)
    cards = [
        (summary_grid[0, 0:2], "Reads", stats["n"], "mapped reads", 0),
        (summary_grid[0, 2:4], "Mean length", stats["mean"], "bp", 1),
        (summary_grid[0, 4:6], "Median length", stats["median"], "bp", 1),
        (summary_grid[1, 1:3], "N50", stats["n50"], "bp", 0),
        (summary_grid[1, 3:5], "Total yield", stats["total_bases"], "bp", 0),
    ]
    for grid_slice, label, value, unit, precision in cards:
        card_ax = fig.add_subplot(grid_slice)
        add_metric_card(card_ax, label, value, unit, precision)

    fig.suptitle("Read QC", fontsize=19, fontweight="bold", color=NAVY, y=0.93)
    fig.savefig(output_file, dpi=300, bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


def main():
    args = parse_args()
    bam_file = args.bam.expanduser()
    output_dir = args.outdir.expanduser()
    output_file = output_dir / args.filename

    if not bam_file.is_file():
        raise FileNotFoundError(f"BAM file does not exist: {bam_file}")

    output_dir.mkdir(parents=True, exist_ok=True)
    lengths = extract_read_lengths(bam_file)
    stats = calculate_statistics(lengths)
    make_figure(lengths, stats, output_file)

    print(f"Input BAM: {bam_file.resolve()}")
    print(f"Mapped reads: {stats['n']:,}")
    print(f"Mean read length: {stats['mean']:.1f} bp")
    print(f"Median read length: {stats['median']:.1f} bp")
    print(f"N50: {stats['n50']:,} bp")
    print(f"Total yield: {stats['total_bases']:,} bp")
    print(f"Saved figure: {output_file.resolve()}")


if __name__ == "__main__":
    main()
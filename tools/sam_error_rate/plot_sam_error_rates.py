
#!/usr/bin/env python3
"""Generate plots from sam_error_rates.py output."""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def normalize_line(line):
    """Normalize whitespace so tabs and spaces are interchangeable."""
    return " ".join(line.split())


def remove_top_right_spines(ax):
    """Remove the top and right borders from a plot."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    parser = argparse.ArgumentParser(
        description="Generate plots from SAM alignment error analysis."
    )
    parser.add_argument("sam", help="Input SAM file")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument(
        "--outdir",
        default="outputs/sam_error_rate",
        help="Directory for generated plots",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    analyzer = script_dir / "sam_error_rates.py"

    # Run the existing analysis and capture its report.
    result = subprocess.run(
        [
            sys.executable,
            str(analyzer),
            args.sam,
            "--reference",
            args.reference,
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    report = result.stdout
    lines = report.splitlines()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def find_heading(heading, start=0):
        """Find a heading by normalized prefix; return its line index."""
        target = normalize_line(heading)

        for i in range(start, len(lines)):
            if normalize_line(lines[i]).startswith(target):
                return i

        raise ValueError(f"Could not find report heading: {heading!r}")

    def section(start, end):
        """Extract lines between two report headings."""
        start_idx = find_heading(start) + 1
        end_idx = find_heading(end, start_idx)

        if end_idx < start_idx:
            raise ValueError(
                f"Invalid report section: {start!r} to {end!r}"
            )

        return lines[start_idx:end_idx]

    def parse_table(table_lines, columns=2):
        """Parse whitespace-delimited tables, skipping headers and blanks."""
        rows = []

        for line in table_lines:
            fields = line.split()

            if len(fields) < columns:
                continue

            try:
                count = int(fields[1])
            except ValueError:
                continue

            rows.append((fields[0], count))

        return rows

    # 1. Substitution spectrum
    sub_rows = parse_table(
        section(
            "Type Count Fraction of substitutions",
            "Transitions total:",
        )
    )

    if not sub_rows:
        raise ValueError("No substitution rows were found in the report.")

    labels = [row[0] for row in sub_rows]
    counts = [row[1] for row in sub_rows]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(labels, counts)
    ax.set(
        title="Substitution spectrum",
        xlabel="Reference base > read base",
        ylabel="Number of substitutions",
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(axis="x", rotation=45)
    remove_top_right_spines(ax)
    fig.tight_layout()
    fig.savefig(outdir / "substitution_spectrum.png", dpi=200)
    plt.close(fig)

    # 2. Indel size distribution (insertions and deletions combined)
    indel_rows = parse_table(
        section(
            "Indel size (bp) Events",
            "Inserted nucleotide Count",
        )
    )

    if indel_rows:
        sizes = [int(row[0]) for row in indel_rows]
        events = [row[1] for row in indel_rows]

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.bar([str(size) for size in sizes], events)
        ax.set(
            title="Indel size distribution",
            xlabel="Indel size (bp)",
            ylabel="Number of events",
        )
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        remove_top_right_spines(ax)
        fig.tight_layout()
        fig.savefig(outdir / "indel_size_distribution.png", dpi=200)
        plt.close(fig)
    else:
        print(
            "Warning: no indel size rows found; "
            "skipping indel size plot."
        )

    # 3. Inserted versus deleted nucleotide composition
    inserted = parse_table(
        section(
            "Inserted nucleotide Count",
            "Deleted reference nucleotide Count",
        )
    )

    deleted = parse_table(
        section(
            "Deleted reference nucleotide Count",
            "Interpretation note:",
        )
    )

    bases = ["A", "C", "G", "T"]
    ins = dict(inserted)
    dels = dict(deleted)

    x = list(range(len(bases)))
    width = 0.36

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(
        [i - width / 2 for i in x],
        [ins.get(base, 0) for base in bases],
        width,
        label="Inserted bases",
    )
    ax.bar(
        [i + width / 2 for i in x],
        [dels.get(base, 0) for base in bases],
        width,
        label="Deleted reference bases",
    )
    ax.set_xticks(x, bases)
    ax.set(
        title="Indel nucleotide composition",
        xlabel="Nucleotide",
        ylabel="Count",
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend()
    remove_top_right_spines(ax)
    fig.tight_layout()
    fig.savefig(outdir / "indel_base_composition.png", dpi=200)
    plt.close(fig)

    # 4. Summary error rates
    rate_names = [
        ("Substitution rate", "Substitutions"),
        ("Indel-event rate", "Indel events"),
        ("Indel-base rate", "Indel bases"),
    ]

    rate_labels = []
    rates = []

    for report_name, plot_label in rate_names:
        match = re.search(
            rf"^{re.escape(report_name)}:\s*"
            rf"[\d,]+\s*/\s*[\d,]+\s*=\s*([\d.]+)%",
            report,
            flags=re.MULTILINE,
        )

        if match is None:
            raise ValueError(
                f"Could not parse {report_name!r} from analyzer output. "
                "Check the report format before generating the summary plot."
            )

        rate_labels.append(plot_label)
        rates.append(float(match.group(1)))

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(rate_labels, rates)

    ax.bar_label(bars, fmt="%.3f%%", padding=3)
    ax.set(
        title="Estimated error rates",
        ylabel="Rate (%)",
    )
    ax.margins(y=0.15)
    remove_top_right_spines(ax)
    fig.tight_layout()
    fig.savefig(outdir / "error_rate_summary.png", dpi=200)
    plt.close(fig)

    print(f"Plots saved to: {outdir.resolve()}")
    for path in sorted(outdir.glob("*.png")):
        print(f"  {path.name}")


if __name__ == "__main__":
    main()

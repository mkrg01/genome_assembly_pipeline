#!/usr/bin/env python3

import argparse
import gzip
import os
from pathlib import Path

cache_root = Path("/tmp") / f"plot_self_dotplot_cache_{os.getuid()}"
cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_root / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter


MIN_ALIGNMENT_LENGTH = 10_000
MIN_IDENTITY = 90.0
DIAGONAL_TOLERANCE_BP = 1_000
ALIGNMENT_COLOR = "#111111"
BOUNDARY_COLOR = "#A8A8A8"
BOUNDARY_LINESTYLE = "dotted"


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def parse_fasta_lengths(fasta_path):
    contigs = []
    seen = set()
    current_name = None
    current_length = 0

    with open_text(fasta_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    contigs.append((current_name, current_length))
                current_name = line[1:].strip().split()[0]
                if not current_name:
                    raise ValueError(f"Empty FASTA record name in {fasta_path}")
                if current_name in seen:
                    raise ValueError(f"Duplicate FASTA record name: {current_name}")
                seen.add(current_name)
                current_length = 0
            else:
                if current_name is None:
                    raise ValueError(f"Sequence data found before a FASTA header in {fasta_path}")
                current_length += len(line.strip())

    if current_name is not None:
        contigs.append((current_name, current_length))

    return contigs


def build_contig_offsets(contigs):
    offsets = {}
    rows = []
    offset = 0
    for name, length in contigs:
        end = offset + length
        offsets[name] = offset
        rows.append((name, length, offset, end))
        offset = end
    return offsets, rows, offset


def write_contig_table(contig_rows, output_path):
    with output_path.open("w") as handle:
        handle.write("name\tlength\tstart\tend\n")
        for name, length, start, end in contig_rows:
            handle.write(f"{name}\t{length}\t{start}\t{end}\n")


def is_main_diagonal(qname, tname, qstart, qend, tstart, tend, strand):
    return (
        strand == "+"
        and qname == tname
        and abs(qstart - tstart) <= DIAGONAL_TOLERANCE_BP
        and abs(qend - tend) <= DIAGONAL_TOLERANCE_BP
    )


def parse_paf(paf_path, offsets):
    forward_segments = []
    reverse_segments = []
    diagonal_segments = []
    total_records = 0
    kept_records = 0

    with paf_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(f"Invalid PAF record at {paf_path}:{line_number}")

            total_records += 1
            qname = fields[0]
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            tname = fields[5]
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            alignment_length = int(fields[10])

            if alignment_length < MIN_ALIGNMENT_LENGTH:
                continue
            if alignment_length == 0:
                continue
            identity = matches / alignment_length * 100
            if identity < MIN_IDENTITY:
                continue
            if qname not in offsets or tname not in offsets:
                raise ValueError(
                    f"PAF record references a sequence not found in the FASTA: "
                    f"{qname} vs {tname}"
                )

            q_abs_start = offsets[qname] + qstart
            q_abs_end = offsets[qname] + qend
            t_abs_start = offsets[tname] + tstart
            t_abs_end = offsets[tname] + tend

            if strand == "+":
                segment = [(t_abs_start, q_abs_start), (t_abs_end, q_abs_end)]
            elif strand == "-":
                segment = [(t_abs_end, q_abs_start), (t_abs_start, q_abs_end)]
            else:
                raise ValueError(f"Invalid strand in {paf_path}:{line_number}: {strand}")

            if is_main_diagonal(qname, tname, qstart, qend, tstart, tend, strand):
                diagonal_segments.append(segment)
            elif strand == "+":
                forward_segments.append(segment)
            else:
                reverse_segments.append(segment)
            kept_records += 1

    return {
        "forward": forward_segments,
        "reverse": reverse_segments,
        "diagonal": diagonal_segments,
        "total_records": total_records,
        "kept_records": kept_records,
    }


def add_line_collection(ax, segments, color, linewidth, alpha, zorder, linestyle="solid"):
    if not segments:
        return
    collection = LineCollection(
        segments,
        colors=color,
        linewidths=linewidth,
        linestyles=linestyle,
        alpha=alpha,
        rasterized=True,
        zorder=zorder,
    )
    ax.add_collection(collection)


def build_boundary_segments(contig_rows, total_length):
    segments = []
    for _name, _length, start, end in contig_rows:
        if start > 0:
            segments.append([(start, 0), (start, total_length)])
            segments.append([(0, start), (total_length, start)])
    return segments


def format_mb(value, _position):
    return f"{value / 1_000_000:g}"


def figure_size(contig_count):
    return min(20.0, max(8.0, 7.0 + contig_count * 0.12))


def contig_label_fontsize(contig_count):
    return max(4.0, min(8.0, 10.0 - contig_count * 0.1))


def add_contig_label_axes(ax, contig_rows):
    if not contig_rows:
        return

    midpoints = [(start + end) / 2 for name, length, start, end in contig_rows]
    labels = [name for name, length, start, end in contig_rows]
    fontsize = contig_label_fontsize(len(contig_rows))

    top_axis = ax.secondary_xaxis("top")
    top_axis.set_xticks(midpoints)
    top_axis.set_xticklabels(
        labels,
        rotation=90,
        fontsize=fontsize,
        ha="left",
        va="center",
        rotation_mode="anchor",
    )
    top_axis.tick_params(length=2, pad=3)
    top_axis.set_xlabel("")

    right_axis = ax.secondary_yaxis("right")
    right_axis.set_yticks(midpoints)
    right_axis.set_yticklabels(labels, fontsize=fontsize)
    right_axis.tick_params(length=2, pad=3)
    right_axis.set_ylabel("")


def escape_mathtext(text):
    replacements = {
        "\\": r"\backslash{}",
        "{": r"\{",
        "}": r"\}",
        "$": r"\$",
        "#": r"\#",
        "%": r"\%",
        "&": r"\&",
        "_": r"\_",
    }
    return "".join(replacements.get(character, character) for character in text)


def format_length_label(length):
    if length >= 1_000_000:
        return f"{length / 1_000_000:g} Mb"
    if length >= 1_000:
        return f"{length / 1_000:g} kb"
    return f"{length:,} bp"


def genome_axis_label(organism_name, min_long_contig_length=None):
    unit_label = "Mb"
    if min_long_contig_length is not None:
        unit_label = f"Mb; contigs >= {format_length_label(min_long_contig_length)}"

    display_name = " ".join(organism_name.replace("_", " ").split())
    if not display_name:
        return f"Genome ({unit_label})"
    italic_name = r"\ ".join(escape_mathtext(word) for word in display_name.split())
    return rf"$\it{{{italic_name}}}$ genome ({unit_label})"


def plot_dotplot(
    contig_rows,
    total_length,
    paf_segments,
    output_pdf,
    axis_label,
):
    size = figure_size(len(contig_rows))
    fig, ax = plt.subplots(figsize=(size, size))

    if total_length == 0:
        ax.text(0.5, 0.5, "No long contigs", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    boundary_segments = build_boundary_segments(contig_rows, total_length)
    add_line_collection(ax, boundary_segments, BOUNDARY_COLOR, 0.35, 0.82, 1, BOUNDARY_LINESTYLE)
    add_line_collection(ax, paf_segments["diagonal"], ALIGNMENT_COLOR, 0.50, 0.70, 2)
    add_line_collection(ax, paf_segments["forward"], ALIGNMENT_COLOR, 0.45, 0.58, 3)
    add_line_collection(ax, paf_segments["reverse"], ALIGNMENT_COLOR, 0.45, 0.58, 3)

    if paf_segments["kept_records"] == 0:
        ax.text(
            0.5,
            0.5,
            "No alignments passed filters",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="#555555",
        )

    ax.set_xlim(0, total_length)
    ax.set_ylim(0, total_length)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(axis_label)
    ax.set_ylabel(axis_label)
    ax.xaxis.set_major_formatter(FuncFormatter(format_mb))
    ax.yaxis.set_major_formatter(FuncFormatter(format_mb))
    add_contig_label_axes(ax, contig_rows)
    ax.tick_params(length=3, width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    fig.tight_layout()
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Plot a minimap2 PAF self-alignment dot plot.")
    parser.add_argument("--fasta", type=Path, required=True, help="Long-contig FASTA.")
    parser.add_argument("--paf", type=Path, required=True, help="Self-alignment PAF from minimap2.")
    parser.add_argument("--pdf", type=Path, required=True, help="Output PDF path.")
    parser.add_argument("--contigs", type=Path, required=True, help="Output contig coordinate TSV path.")
    parser.add_argument("--organism-name", required=True, help="Organism name used for axis labels.")
    parser.add_argument(
        "--min-long-contig-length",
        type=int,
        default=None,
        help="Minimum contig length included in the plotted FASTA.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    contigs = parse_fasta_lengths(args.fasta)
    offsets, contig_rows, total_length = build_contig_offsets(contigs)
    write_contig_table(contig_rows, args.contigs)
    paf_segments = parse_paf(args.paf, offsets)
    plot_dotplot(
        contig_rows,
        total_length,
        paf_segments,
        args.pdf,
        genome_axis_label(args.organism_name, args.min_long_contig_length),
    )
    print(
        f"Plotted {paf_segments['kept_records']} of {paf_segments['total_records']} "
        f"PAF records with length >= {MIN_ALIGNMENT_LENGTH} bp and identity >= {MIN_IDENTITY:.1f}%."
    )


if __name__ == "__main__":
    main()

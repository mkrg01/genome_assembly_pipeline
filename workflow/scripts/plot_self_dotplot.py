#!/usr/bin/env python3

import argparse
import gzip
import json
import os
from pathlib import Path

cache_root = Path("/tmp") / f"plot_self_dotplot_cache_{os.getuid()}"
cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_root / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter

from plot_style import apply_matplotlib_style

apply_matplotlib_style()


MIN_ALIGNMENT_LENGTH = 100_000
ALIGNMENT_COLOR = "#222222"
BOUNDARY_COLOR = "#CBCBCB"
BOUNDARY_LINESTYLE = "dotted"
CONTIG_LABEL_SIZE = 7.0
LABEL_GAP_PT = 0.75
LABEL_LEADER_PT = 12.0
LABEL_COLUMN_GAP_PT = 4.0


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
        and qstart == tstart
        and qend == tend
    )


def parse_paf(paf_path, offsets, min_alignment_length=MIN_ALIGNMENT_LENGTH):
    segments = []
    total_records = 0
    counts = dict.fromkeys(("short_span", "self_diagonal"), 0)

    with open_text(paf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            # Keep optional tags together: CIGAR strings can be very large.
            fields = line.split("\t", 12)
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
            if qname not in offsets or tname not in offsets:
                raise ValueError(
                    f"PAF record references a sequence not found in the FASTA: "
                    f"{qname} vs {tname}"
                )
            if strand not in ("+", "-"):
                raise ValueError(f"Invalid strand in {paf_path}:{line_number}: {strand}")
            if not (0 <= qstart < qend <= int(fields[1]) and 0 <= tstart < tend <= int(fields[6])):
                raise ValueError(f"Invalid coordinates in {paf_path}:{line_number}")
            # Use coordinate spans only. Mapping-only PAF block lengths and
            # matching-base counts are approximate and do not define identity.
            if min(qend - qstart, tend - tstart) < min_alignment_length:
                counts["short_span"] += 1
                continue
            if is_main_diagonal(qname, tname, qstart, qend, tstart, tend, strand):
                counts["self_diagonal"] += 1
                continue
            q_abs_start = offsets[qname] + qstart
            q_abs_end = offsets[qname] + qend
            t_abs_start = offsets[tname] + tstart
            t_abs_end = offsets[tname] + tend

            if strand == "+":
                segment = [(t_abs_start, q_abs_start), (t_abs_end, q_abs_end)]
            else:
                segment = [(t_abs_end, q_abs_start), (t_abs_start, q_abs_end)]

            segments.append(segment)

    return {
        "segments": segments,
        "total_records": total_records,
        "kept_records": len(segments),
        "filtered_records": counts,
    }


def add_line_collection(ax, segments, color, linewidth, alpha, zorder, linestyle="solid", capstyle="butt"):
    if not segments:
        return
    collection = LineCollection(
        segments,
        colors=color,
        linewidths=linewidth,
        linestyles=linestyle,
        capstyle=capstyle,
        alpha=alpha,
        rasterized=False,
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
    # Keep a fixed publication width so additional contigs do not shrink text.
    return 180 / 25.4


def contig_label_fontsize(contig_count):
    return CONTIG_LABEL_SIZE


def pack_label_positions(desired, spacing, lower, upper):
    """Least-squares ordered placement with prescribed gaps and bounded centers."""
    if not desired:
        return []
    gaps = [spacing] * (len(desired) - 1) if isinstance(spacing, (int, float)) else list(spacing)
    if len(gaps) != len(desired) - 1 or any(gap <= 0 for gap in gaps) or any(a > b for a, b in zip(desired, desired[1:])):
        raise ValueError("Label positions must be ordered and spacing must be positive")
    offsets = [0]
    for gap in gaps:
        offsets.append(offsets[-1] + gap)
    ceiling = upper - offsets[-1]
    if ceiling < lower - 1e-9:
        raise ValueError("Insufficient space for contig labels")
    ceiling = max(lower, ceiling)
    # Subtract the required gaps, then use pooled adjacent violators to find
    # the closest nondecreasing sequence. No random repulsion or name thinning.
    blocks = []
    for index, position in enumerate(desired):
        blocks.append((position - offsets[index], 1))
        while len(blocks) > 1:
            left, right = blocks[-2:]
            if left[0] / left[1] <= right[0] / right[1]:
                break
            blocks[-2:] = [(left[0] + right[0], left[1] + right[1])]
    fitted = []
    for total, count in blocks:
        fitted.extend([min(ceiling, max(lower, total / count))] * count)
    return [position + offset for position, offset in zip(fitted, offsets)]


def add_contig_label_axes(ax, contig_rows):
    """Place literal names outside the plot, with ordered leaders when needed."""
    if not contig_rows:
        return
    fig = ax.figure
    labels = [name for name, length, start, end in contig_rows]
    fontsize = contig_label_fontsize(len(contig_rows))
    fig.tight_layout(pad=0.8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    points_per_pixel = 72 / fig.dpi
    probes = [fig.text(0, 0, name, fontsize=fontsize, parse_math=False) for name in labels]
    boxes = [text.get_window_extent(renderer) for text in probes]
    text_width = max(box.width for box in boxes) * points_per_pixel
    text_height = max(box.height for box in boxes) * points_per_pixel
    # PDF text can occupy the full font ascent/descent beyond Agg's glyph box.
    font = font_manager.get_font(font_manager.findfont(probes[0].get_fontproperties()))
    text_height = max(text_height, fontsize * (font.ascender - font.descender) / font.units_per_EM + 0.1)
    for text in probes:
        text.remove()

    width, height = fig.get_size_inches() * 72
    left, bottom = ax.bbox.x0 * points_per_pixel, ax.bbox.y0 * points_per_pixel
    total_length = contig_rows[-1][3]
    midpoints = [(start + end) / (2 * total_length) for _, _, start, end in contig_rows]
    candidates = []
    for lanes in (1, 2):
        margin = LABEL_LEADER_PT + 3 + lanes * text_width + (lanes - 1) * LABEL_COLUMN_GAP_PT + 4
        side = min(width - left - margin, height - bottom - margin)
        if side < 100 / 25.4 * 72:
            continue
        desired = [midpoint * side for midpoint in midpoints]
        pitch = text_height + LABEL_GAP_PT
        lane_indices = [0] * len(labels)
        if lanes == 2:
            crowded = [False] * len(labels)
            for index in range(1, len(labels)):
                if desired[index] - desired[index - 1] < pitch:
                    crowded[index - 1] = crowded[index] = True
            run = 0
            for index, collision in enumerate(crowded):
                lane_indices[index] = run % 2 if collision else 0
                run = run + 1 if collision else 0
        gaps = [pitch if a == b else pitch / 2 for a, b in zip(lane_indices, lane_indices[1:])]
        try:
            positions = pack_label_positions(
                desired, gaps,
                text_height / 2,
                side - text_height / 2,
            )
        except ValueError:
            continue
        shift = max(abs(position - target) for position, target in zip(positions, desired))
        candidates.append((shift, side, positions, lane_indices))
        if shift <= 30:
            break
    if not candidates:
        raise ValueError(
            "All contig labels cannot fit at 7 pt in the fixed-width figure. "
            "Increase the figure dimensions; no names have been omitted."
        )
    _, side, positions, lane_indices = min(candidates, key=lambda candidate: candidate[0])
    ax.set_position([left / width, bottom / height, side / width, side / height])

    ticks, leaders = [], []
    for index, (name, midpoint, position) in enumerate(zip(labels, midpoints, positions)):
        # Interleaved ports keep outer leaders in the gaps between inner labels.
        lane = lane_indices[index]
        center = position / side
        offset = LABEL_LEADER_PT + 3 + lane * (text_width + LABEL_COLUMN_GAP_PT)
        for edge in ("top", "right"):
            if edge == "top":
                xy = (center, 1 + offset / side)
                tick = [(midpoint, 1), (midpoint, 1 + 2 / side)]
                leader = [tick[-1], (center, 1 + (LABEL_LEADER_PT - 2) / side),
                          (center, 1 + (offset - 1.5) / side)]
            else:
                xy = (1 + offset / side, center)
                tick = [(1, midpoint), (1 + 2 / side, midpoint)]
                leader = [tick[-1], (1 + (LABEL_LEADER_PT - 2) / side, center),
                          (1 + (offset - 1.5) / side, center)]
            text = ax.text(
                *xy, name, transform=ax.transAxes, fontsize=fontsize,
                ha="left", va="center", rotation=90 if edge == "top" else 0,
                rotation_mode="anchor", clip_on=False, parse_math=False,
            )
            text.set_gid(f"contig-label-{edge}-{index}")
            ticks.append(tick)
            if lane or abs(position - midpoint * side) > 0.25:
                leaders.append(leader)
    for segments, color, linewidth, gid in (
        (ticks, "black", 0.6, "contig-label-ticks"),
        (leaders, "#777777", 0.35, "contig-label-leaders"),
    ):
        collection = LineCollection(
            segments, colors=color, linewidths=linewidth, transform=ax.transAxes,
            clip_on=False, zorder=4,
        )
        collection.set_gid(gid)
        ax.add_collection(collection)


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
    min_alignment_length=MIN_ALIGNMENT_LENGTH,
):
    size = figure_size(len(contig_rows))
    fig, ax = plt.subplots(figsize=(size, size), dpi=144)

    if total_length == 0:
        ax.text(0.5, 0.5, "No long contigs", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        fig.savefig(output_pdf, dpi=600)
        plt.close(fig)
        return

    boundary_segments = build_boundary_segments(contig_rows, total_length)
    add_line_collection(ax, boundary_segments, BOUNDARY_COLOR, 0.25, 0.65, 1, BOUNDARY_LINESTYLE)
    # Mark identical coordinates as a guide, separate from homology matches.
    ax.plot(
        [0, total_length], [0, total_length], color=BOUNDARY_COLOR,
        linestyle="solid", linewidth=0.5, zorder=1,
    )
    # Flat caps make a sublinewidth match look perpendicular to its true
    # direction. Round caps represent unresolved matches as dots instead.
    add_line_collection(ax, paf_segments["segments"], ALIGNMENT_COLOR, 0.4, 1.0, 3, capstyle="round")

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
    ax.set_xlabel(axis_label, fontsize=9)
    ax.set_ylabel(axis_label, fontsize=9)
    ax.xaxis.set_major_formatter(FuncFormatter(format_mb))
    ax.yaxis.set_major_formatter(FuncFormatter(format_mb))
    ax.tick_params(length=3, width=0.8, labelsize=7)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    try:
        add_contig_label_axes(ax, contig_rows)
        fig.savefig(output_pdf, dpi=600)
    finally:
        plt.close(fig)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Plot a minimap2 PAF self-alignment dot plot.")
    parser.add_argument("--fasta", type=Path, required=True, help="Long-contig FASTA.")
    parser.add_argument("--paf", type=Path, required=True, help="Self-alignment PAF from minimap2.")
    parser.add_argument("--pdf", type=Path, required=True, help="Output PDF path.")
    parser.add_argument("--contigs", type=Path, required=True, help="Output contig coordinate TSV path.")
    parser.add_argument("--summary", type=Path, help="Optional output JSON with filtering counts.")
    parser.add_argument("--organism-name", required=True, help="Organism name used for axis labels.")
    parser.add_argument(
        "--min-long-contig-length",
        type=int,
        default=None,
        help="Minimum contig length included in the plotted FASTA.",
    )
    parser.add_argument(
        "--min-alignment-length", type=int, default=MIN_ALIGNMENT_LENGTH,
        help=f"Minimum coordinate span on each sequence in bp (default: {MIN_ALIGNMENT_LENGTH}).",
    )
    args = parser.parse_args(argv)
    if args.min_alignment_length < 1:
        parser.error("--min-alignment-length must be a positive integer")
    return args


def main():
    args = parse_args()
    contigs = parse_fasta_lengths(args.fasta)
    offsets, contig_rows, total_length = build_contig_offsets(contigs)
    write_contig_table(contig_rows, args.contigs)
    paf_segments = parse_paf(args.paf, offsets, args.min_alignment_length)
    plot_dotplot(
        contig_rows,
        total_length,
        paf_segments,
        args.pdf,
        genome_axis_label(args.organism_name, args.min_long_contig_length),
        min_alignment_length=args.min_alignment_length,
    )
    summary = {
        "total_records": paf_segments["total_records"],
        "plotted_records": paf_segments["kept_records"],
        "filtered_records": paf_segments["filtered_records"],
        "min_alignment_length": args.min_alignment_length,
        "contigs": len(contigs),
        "total_length": total_length,
    }
    if args.summary:
        args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"Plotted {paf_segments['kept_records']} of {paf_segments['total_records']} "
        f"PAF records with spans >= {args.min_alignment_length} bp on each sequence."
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()

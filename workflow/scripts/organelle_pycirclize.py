import argparse
import math
import os
import tempfile
from collections import Counter
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

FEATURE_ORDER = ("CDS", "rRNA", "tRNA")
REPEAT_FEATURE_TYPE = "repeat_region"
FEATURE_COLORS = {
    "gene": "#4C72B0",
    "repeat": "#B7B7B7",
}
OUTER_STRAND = 1
INNER_STRAND = -1
SECTOR_CLOCKWISE = False
BASELINE_R = 78.0
FEATURE_WIDTH = 4.0
OUTER_FEATURE_R_LIM = (BASELINE_R, BASELINE_R + FEATURE_WIDTH)
INNER_FEATURE_R_LIM = (BASELINE_R - FEATURE_WIDTH, BASELINE_R)
LABEL_LINE_GAP = 0.4
OUTER_LABEL_R = OUTER_FEATURE_R_LIM[1] + 2.2
INNER_LABEL_R = INNER_FEATURE_R_LIM[0] - 3.0
TITLE_LINE_AXES_OFFSET = 0.035
DIRECTION_ARROW_OUTER_R_LIM = (86.2, 87.4)
DIRECTION_ARROW_INNER_R_LIM = (68.2, 69.4)
DIRECTION_ARROW_LENGTH_FRACTION = 0.018
DIRECTION_ARROW_MAX_LENGTH = 3_000
DIRECTION_ARROW_MIN_LENGTH = 600
DIRECTION_ARROW_COLOR = "#B8B8B8"
DIRECTION_ARROW_ALPHA = 0.5
DIRECTION_ARROW_DISPLAY_DEG = 5.0
REPEAT_BASELINE_R = 49.0
REPEAT_R_LIM = (48.2, 49.8)
REPEAT_LABEL_R = 51.6


def first_qualifier(feature, keys):
    qualifiers = getattr(feature, "qualifiers", {})
    for key in keys:
        values = qualifiers.get(key)
        if values:
            return str(values[0])
    return None


def short_label(label, max_len=28):
    if label is None:
        return None
    label = " ".join(str(label).split())
    if len(label) <= max_len:
        return label
    return f"{label[: max_len - 1]}."


def feature_label(feature):
    label = first_qualifier(feature, ("gene", "locus_tag", "product", "note"))
    return short_label(label)


def repeat_label(feature):
    label = first_qualifier(feature, ("note", "rpt_type"))
    if label is None:
        return "repeat"

    words = label.split()
    if len(words) >= 3 and " ".join(words[:2]).lower() == "inverted repeat":
        repeat_name = words[-1]
        if len(repeat_name) == 1 and repeat_name.isalpha():
            return f"IR{repeat_name.lower()}"
    return short_label(label, max_len=18)


def feature_parts(feature):
    parts = getattr(feature.location, "parts", None)
    if parts is None:
        parts = [feature.location]

    segments = []
    for part in parts:
        start = int(part.start)
        end = int(part.end)
        strand = part.strand if part.strand is not None else feature.location.strand
        if strand is None:
            strand = INNER_STRAND
        segments.append((start, end, strand))
    return segments


def feature_part_labels(feature):
    label = feature_label(feature)
    if label is None:
        return []

    parts = feature_parts(feature)
    part_count = len(parts)
    labels = []
    for index, (start, end, strand) in enumerate(parts, start=1):
        part_label = label
        if part_count > 1:
            part_label = f"{label} exon {index}/{part_count}"
        labels.append(((start + end) / 2, strand, part_label))
    return labels


def longest_part_anchor(feature):
    start, end, strand = max(feature_parts(feature), key=lambda part: part[1] - part[0])
    return (start + end) / 2, strand


def feature_r_lim(strand):
    if strand == OUTER_STRAND:
        return OUTER_FEATURE_R_LIM, FEATURE_COLORS["gene"]
    return INNER_FEATURE_R_LIM, FEATURE_COLORS["gene"]


def feature_label_layout(strand):
    if strand == OUTER_STRAND:
        return {
            "line_r_lim": (OUTER_FEATURE_R_LIM[1], OUTER_LABEL_R - LABEL_LINE_GAP),
            "label_r": OUTER_LABEL_R,
        }
    return {
        "line_r_lim": (INNER_LABEL_R + LABEL_LINE_GAP, INNER_FEATURE_R_LIM[0]),
        "label_r": INNER_LABEL_R,
    }


def is_right_rad(rad):
    deg = math.degrees(rad)
    return -360 <= deg < -180 or 0 <= deg < 180


def feature_label_text_params(rad, strand):
    is_right = is_right_rad(rad)
    outer = strand == OUTER_STRAND
    if outer:
        ha = "left" if is_right else "right"
    else:
        ha = "right" if is_right else "left"
    return {
        "ha": ha,
        "va": "center_baseline",
        "rotation": 90 - math.degrees(rad) if is_right else 270 - math.degrees(rad),
        "rotation_mode": "anchor",
    }


def display_assembly_name(assembly_name):
    return str(assembly_name).replace("_", " ")


def genome_size_label(seqid2size):
    total_size = sum(seqid2size.values())
    return f"{total_size:,} bp"


def organelle_genome_label(organelle):
    organelle_labels = {
        "chloroplast": "chloroplast genome",
        "mitochondrion": "mitochondrial genome",
    }
    if organelle is None:
        return "genome"
    return organelle_labels.get(organelle, f"{organelle} genome")


def map_title(assembly_name, organelle, seqid2size=None):
    title = f"{display_assembly_name(assembly_name)}\n{organelle_genome_label(organelle)}"
    if seqid2size is not None:
        title = f"{title}\n{genome_size_label(seqid2size)}"
    return title


def map_title_lines(assembly_name, organelle, seqid2size):
    return (
        display_assembly_name(assembly_name),
        organelle_genome_label(organelle),
        genome_size_label(seqid2size),
    )


def add_map_title(circos, assembly_name, organelle, seqid2size, title=None):
    if title:
        circos.text(title, size=11)
        return

    species_name, genome_label, size_label = map_title_lines(
        assembly_name, organelle, seqid2size
    )

    def plot_title(ax):
        common_kwargs = {
            "transform": ax.transAxes,
            "ha": "center",
            "va": "center",
            "size": 11,
        }
        ax.text(
            0.5,
            0.5 + TITLE_LINE_AXES_OFFSET,
            species_name,
            fontstyle="italic",
            **common_kwargs,
        )
        ax.text(0.5, 0.5, genome_label, **common_kwargs)
        ax.text(0.5, 0.5 - TITLE_LINE_AXES_OFFSET, size_label, **common_kwargs)

    circos._plot_funcs.append(plot_title)


def load_genbank_records(genbank_path):
    from Bio import SeqIO

    records = list(SeqIO.parse(str(genbank_path), "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found in {genbank_path}")
    return records


def record_seqid(record):
    for value in (getattr(record, "id", None), getattr(record, "name", None)):
        if value and value != "<unknown id>":
            return str(value)
    return str(record.description)


def records_to_seqid2size(records):
    return {record_seqid(record): len(record.seq) for record in records}


def get_seqid2features(records, feature_type):
    seqid2features = {}
    for record in records:
        features = [
            feature
            for feature in record.features
            if feature.type == feature_type and feature.location is not None
        ]
        seqid2features[record_seqid(record)] = features
    return seqid2features


def direction_arrow_length(span):
    return min(
        max(span * DIRECTION_ARROW_LENGTH_FRACTION, DIRECTION_ARROW_MIN_LENGTH),
        DIRECTION_ARROW_MAX_LENGTH,
    )


def sector_x_at_display_deg(sector, display_deg):
    if sector is None:
        return None

    target_rad = math.radians(display_deg) % (2 * math.pi)
    rad_start, rad_end = min(sector.rad_lim), max(sector.rad_lim)
    if not rad_start <= target_rad <= rad_end:
        return None

    size_ratio = sector.rad_size / sector.size if sector.size else 0
    if size_ratio == 0:
        return None

    x = sector.start + ((target_rad - rad_start) / size_ratio)
    if not sector.clockwise:
        x = (sector.start + sector.end) - x
    if not sector.start <= x <= sector.end:
        return None
    return x


def direction_arrow_bounds(track):
    start = track.start
    end = track.end
    span = end - start
    if span <= 0:
        return None

    arrow_len = direction_arrow_length(span)
    sector = getattr(track, "parent_sector", None)
    center = sector_x_at_display_deg(sector, DIRECTION_ARROW_DISPLAY_DEG)
    if center is None:
        pad = max(span * 0.02, 1)
        center = start + pad + arrow_len / 2

    arrow_start = max(start, min(center - arrow_len / 2, end - arrow_len))
    arrow_end = arrow_start + arrow_len
    if arrow_end <= arrow_start:
        return None
    return arrow_start, arrow_end


def add_direction_arrows(track):
    bounds = direction_arrow_bounds(track)
    if bounds is None:
        return
    arrow_start, arrow_end = bounds

    track.arrow(
        arrow_start,
        arrow_end,
        r_lim=DIRECTION_ARROW_OUTER_R_LIM,
        head_length=2.5,
        shaft_ratio=0.45,
        fc=DIRECTION_ARROW_COLOR,
        ec="none",
        alpha=DIRECTION_ARROW_ALPHA,
    )
    track.arrow(
        arrow_end,
        arrow_start,
        r_lim=DIRECTION_ARROW_INNER_R_LIM,
        head_length=2.5,
        shaft_ratio=0.45,
        fc=DIRECTION_ARROW_COLOR,
        ec="none",
        alpha=DIRECTION_ARROW_ALPHA,
    )


def add_feature_label(track, x, strand, label):
    layout = feature_label_layout(strand)
    text_params = feature_label_text_params(track.x_to_rad(x), strand)
    track._simpleline(
        (x, x),
        layout["line_r_lim"],
        ec="#999999",
        lw=0.25,
    )
    track.text(
        label,
        x=x,
        r=layout["label_r"],
        size=4.5,
        adjust_rotation=False,
        color="#111111",
        **text_params,
    )


def plot_feature_boxes(track, feature):
    for start, end, strand in feature_parts(feature):
        r_lim, color = feature_r_lim(strand)
        track.rect(
            start,
            end,
            r_lim=r_lim,
            fc=color,
            ec=color,
            lw=0.15,
            alpha=0.85,
        )


def plot_repeat_region(track, feature):
    for start, end, _strand in feature_parts(feature):
        track.rect(
            start,
            end,
            r_lim=REPEAT_R_LIM,
            fc=FEATURE_COLORS["repeat"],
            ec=FEATURE_COLORS["repeat"],
            lw=0.15,
            alpha=0.75,
        )


def add_repeat_region_label(track, feature):
    x, _strand = longest_part_anchor(feature)
    track.text(
        repeat_label(feature),
        x=x,
        r=REPEAT_LABEL_R,
        size=4.2,
        orientation="vertical",
        color="#555555",
    )


def add_repeat_regions(sector, track, features, *, label_features=True):
    if not features:
        return

    sector.line(r=REPEAT_BASELINE_R, lw=0.45, color="#999999")
    for feature in features:
        plot_repeat_region(track, feature)
        if label_features:
            add_repeat_region_label(track, feature)


def render_organelle_pycirclize(
    genbank_path,
    output_path,
    *,
    assembly_name=None,
    organelle=None,
    title=None,
    label_features=True,
):
    from pycirclize import Circos

    genbank_path = Path(genbank_path)
    output_path = Path(output_path)

    records = load_genbank_records(genbank_path)
    seqid2size = records_to_seqid2size(records)

    sector2clockwise = {seqid: SECTOR_CLOCKWISE for seqid in seqid2size}
    circos = Circos(
        seqid2size,
        space=0 if len(seqid2size) == 1 else 2,
        sector2clockwise=sector2clockwise,
    )
    add_map_title(
        circos,
        assembly_name or genbank_path.stem,
        organelle,
        seqid2size,
        title=title,
    )

    seqid2feature_type2features = {
        feature_type: get_seqid2features(records, feature_type)
        for feature_type in FEATURE_ORDER
    }
    seqid2repeat_features = get_seqid2features(records, REPEAT_FEATURE_TYPE)
    feature_counts = Counter()

    for sector in circos.sectors:
        if len(seqid2size) > 1:
            sector.text(sector.name, r=108, size=7, orientation="vertical")

        sector.line(r=BASELINE_R, lw=0.55, color="#555555")
        track = sector.add_track((46, 94))
        repeat_features = seqid2repeat_features.get(sector.name, [])
        add_repeat_regions(
            sector,
            track,
            repeat_features,
            label_features=label_features,
        )
        feature_counts[REPEAT_FEATURE_TYPE] += len(repeat_features)
        add_direction_arrows(track)
        seen_labels = set()

        for feature_type in FEATURE_ORDER:
            features = seqid2feature_type2features[feature_type].get(sector.name, [])
            for feature in features:
                plot_feature_boxes(track, feature)
                feature_counts[feature.type] += 1

                if label_features and feature.type in FEATURE_ORDER:
                    for x, strand, label in feature_part_labels(feature):
                        label_key = (round(x, 3), strand, label)
                        if label_key in seen_labels:
                            continue
                        seen_labels.add(label_key)
                        add_feature_label(track, x, strand, label)

    fig = circos.plotfig()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return feature_counts


def parse_args():
    parser = argparse.ArgumentParser(
        description="Draw a pyCirclize circular organelle map from a GenBank file."
    )
    parser.add_argument("--genbank", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--assembly-name")
    parser.add_argument("--organelle")
    parser.add_argument("--title")
    parser.add_argument("--no-labels", action="store_true")
    return parser.parse_args()


def run_from_args(args):
    render_organelle_pycirclize(
        args.genbank,
        args.output,
        assembly_name=args.assembly_name,
        organelle=args.organelle,
        title=args.title,
        label_features=not args.no_labels,
    )


def run_from_snakemake(snakemake):
    label_features = bool(snakemake.params.get("label_features", True))
    render_organelle_pycirclize(
        snakemake.input.annotation,
        snakemake.output[0],
        assembly_name=snakemake.params.get("assembly_name", None),
        organelle=snakemake.params.get("organelle", None),
        title=snakemake.params.get("title", None),
        label_features=label_features,
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
elif __name__ == "__main__":
    run_from_args(parse_args())

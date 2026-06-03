#!/usr/bin/env python3
"""Render a gbdraw circular organelle map from a GenBank file."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import os
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from gbdraw.api import assemble_circular_diagram_from_record, save_figure_to
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.labels.filtering import read_qualifier_priority_file
from gbdraw.render.drawers.circular.definition import DefinitionDrawer

DEFAULT_FEATURES = "CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA"


def parse_formats(value: str) -> list[str]:
    formats: list[str] = []
    for item in value.split(","):
        fmt = item.strip().lower()
        if fmt and fmt not in formats:
            formats.append(fmt)
    if not formats:
        formats.append("pdf")
    return formats


@contextmanager
def circular_accession_visibility(show_accession: bool):
    original_draw = DefinitionDrawer.draw

    if show_accession:
        yield
        return

    def draw_without_accession(
        self,
        definition_group,
        title_x,
        title_y,
        species_parts,
        strain_parts,
        organelle_parts,
        replicon_parts,
        gc_percent,
        accession,
        record_length,
        *args,
        **kwargs,
    ):
        kwargs["show_accession"] = False
        return original_draw(
            self,
            definition_group,
            title_x,
            title_y,
            species_parts,
            strain_parts,
            organelle_parts,
            replicon_parts,
            gc_percent,
            accession,
            record_length,
            *args,
            **kwargs,
        )

    DefinitionDrawer.draw = draw_without_accession
    try:
        yield
    finally:
        DefinitionDrawer.draw = original_draw


def exact_output_paths(output_prefix: Path, formats: Iterable[str]) -> list[Path]:
    return [Path(f"{output_prefix}.{fmt}") for fmt in formats]


def parse_features(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def has_selected_features(record, selected_features: Iterable[str]) -> bool:
    selected = set(selected_features)
    return any(feature.type in selected for feature in record.features)


def save_with_exact_prefix(canvas, formats: list[str], output_prefix: Path, overwrite: bool) -> list[Path]:
    output_dir = output_prefix.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    target_paths = exact_output_paths(output_prefix, formats)
    if not overwrite:
        existing = [str(path) for path in target_paths if path.exists()]
        if existing:
            raise FileExistsError(
                "Output file(s) already exist: " + ", ".join(existing) + ". Use --overwrite."
            )

    safe_base = f"{output_prefix.name.replace('.', '_')}_tmp_gbdraw_{os.getpid()}"
    save_figure_to(
        canvas,
        formats,
        output_dir=str(output_dir),
        output_prefix=safe_base,
        overwrite=True,
    )

    written: list[Path] = []
    for fmt in formats:
        source = output_dir / f"{safe_base}.{fmt}"
        target = Path(f"{output_prefix}.{fmt}")
        if overwrite:
            source.replace(target)
        else:
            source.rename(target)
        written.append(target)
    temporary_svg = output_dir / f"{safe_base}.svg"
    if "svg" not in formats and temporary_svg.exists():
        temporary_svg.unlink()
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render a gbdraw circular organelle map with gene-prioritized labels."
    )
    parser.add_argument("--gbk", required=True, type=Path, help="Input GenBank file.")
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help="Output prefix without the final format extension.",
    )
    parser.add_argument(
        "--formats",
        default="pdf",
        help="Comma-separated output formats.",
    )
    parser.add_argument(
        "--species",
        help="Species name to show in the center. Overrides the GenBank source organism.",
    )
    parser.add_argument(
        "--plain-species",
        action="store_true",
        help="Do not italicize the species name supplied with --species.",
    )
    parser.add_argument(
        "--qualifier-priority",
        type=Path,
        help="Two-column TSV specifying label qualifier priority by feature type.",
    )
    parser.add_argument(
        "--labels",
        choices=("none", "outside", "both"),
        default="both",
        help="Label placement mode.",
    )
    parser.add_argument(
        "--features",
        default=DEFAULT_FEATURES,
        help="Comma-separated GenBank feature types to draw.",
    )
    parser.add_argument(
        "--include-repeat-region",
        action="store_true",
        help="Also draw repeat_region features.",
    )
    parser.add_argument(
        "--single-strand-ring",
        action="store_true",
        help="Draw both strands on one feature ring.",
    )
    parser.add_argument(
        "--track-type",
        choices=("middle", "tuckin", "spreadout"),
        default="middle",
        help="Feature ring placement mode.",
    )
    parser.add_argument(
        "--hide-accession",
        action="store_true",
        help="Hide the record ID/accession in the center definition block.",
    )
    parser.add_argument(
        "--show-gc-plot",
        action="store_true",
        help="Show the GC content track. The center GC summary is kept even when hidden.",
    )
    parser.add_argument(
        "--show-skew",
        action="store_true",
        help="Show the GC skew track. Hidden by default.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    formats = parse_formats(args.formats)
    selected_features = parse_features(args.features)
    if args.include_repeat_region and "repeat_region" not in selected_features:
        selected_features.append("repeat_region")
    record = SeqIO.read(args.gbk, "genbank")

    species = args.species
    if species and not args.plain_species:
        species = f"<i>{species}</i>"

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_gc=bool(args.show_gc_plot),
        show_skew=bool(args.show_skew),
        show_labels=args.labels != "none",
        allow_inner_labels=args.labels == "both",
        strandedness=not args.single_strand_ring,
        track_type=args.track_type,
    )
    if args.qualifier_priority:
        config_dict["labels"]["filtering"]["qualifier_priority_df"] = (
            read_qualifier_priority_file(str(args.qualifier_priority))
        )

    with circular_accession_visibility(not args.hide_accession):
        canvas = assemble_circular_diagram_from_record(
            record,
            config_dict=config_dict,
            selected_features_set=selected_features,
            output_prefix=args.output_prefix.name,
            species=species,
            legend="right" if has_selected_features(record, selected_features) else "none",
        )

    written = save_with_exact_prefix(canvas, formats, args.output_prefix, args.overwrite)
    for path in written:
        print(path)


if __name__ == "__main__":
    main()

import argparse
from pathlib import Path


ORGANELLE_LABELS = {
    "mitochondrion": "mitochondrial",
    "mitochondria": "mitochondrial",
    "mito": "mitochondrial",
    "chloroplast": "chloroplast",
    "plastid": "plastid",
    "pltd": "plastid",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write a README.md describing organelle submission file provenance."
    )
    parser.add_argument("--organelle", required=True)
    parser.add_argument("--annotation-tool")
    parser.add_argument("--assembly-name", required=True)
    parser.add_argument("--assembly-version", required=True)
    parser.add_argument("--input-genome", type=Path, required=True)
    parser.add_argument("--input-annotation", type=Path)
    parser.add_argument("--output-genome", type=Path, required=True)
    parser.add_argument("--output-annotation", type=Path)
    parser.add_argument("--output-readme", type=Path, required=True)
    return parser.parse_args()


def organelle_label(organelle: str):
    try:
        return ORGANELLE_LABELS[organelle]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported organelle '{organelle}'. Expected one of: {', '.join(ORGANELLE_LABELS)}"
        ) from exc


def quote_path(path: Path):
    return f"`{path.as_posix()}`"


def write_readme(args):
    label = organelle_label(args.organelle)
    args.output_readme.parent.mkdir(parents=True, exist_ok=True)
    annotation_tool = getattr(args, "annotation_tool", None)
    input_annotation = getattr(args, "input_annotation", None)
    output_annotation = getattr(args, "output_annotation", None)
    annotation_args = [annotation_tool, input_annotation, output_annotation]
    if any(annotation_args) and not all(annotation_args):
        raise ValueError(
            "--annotation-tool, --input-annotation, and --output-annotation must be "
            "provided together."
        )

    provenance_rows = [
        (
            f"| {quote_path(args.output_genome)} | {quote_path(args.input_genome)} | "
            "The organelle genome FASTA was copied without sequence changes and compressed "
            "with `gzip`. |"
        )
    ]
    notes = [
        "- The organelle genome was staged from the Oatk assembly output selected by "
        "`oatk_organelle` in `config/config.yml`."
    ]
    if annotation_tool and input_annotation and output_annotation:
        provenance_rows.append(
            (
                f"| {quote_path(output_annotation)} | {quote_path(input_annotation)} | "
                f"The organelle annotation file generated with `{annotation_tool}` was copied "
                "without content changes and compressed with `gzip`. |"
            )
        )
        notes.append(
            "- The organelle annotation was staged from the configured organelle annotation "
            "tool selected by `organelle_annotation` in `config/config.yml`."
        )
    else:
        notes.append(
            "- No organelle annotation file was staged because no annotation tool was "
            "configured for this organelle."
        )

    content = f"""# Organelle Submission Files

This directory contains files staged for submission of the {label} genome for `{args.assembly_name}` (version `{args.assembly_version}`).

## File Provenance

| Submission file | Source file | How the submission file was generated |
| --- | --- | --- |
{chr(10).join(provenance_rows)}

## Notes

{chr(10).join(notes)}
- No automatic identifier renaming or annotation rewriting was applied to the organelle inputs.
- This README was generated automatically by `workflow/scripts/write_organelle_submission_readme.py`.
"""

    args.output_readme.write_text(content)


def main():
    args = parse_args()
    write_readme(args)


if __name__ == "__main__":
    main()

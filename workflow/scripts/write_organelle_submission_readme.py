import argparse
from pathlib import Path


ORGANELLE_LABELS = {
    "mito": "mitochondrial",
    "pltd": "plastid",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write a README.md describing organelle submission file provenance."
    )
    parser.add_argument("--organelle", required=True)
    parser.add_argument("--assembly-name", required=True)
    parser.add_argument("--assembly-version", required=True)
    parser.add_argument("--input-genome", type=Path, required=True)
    parser.add_argument("--input-annotation", type=Path, required=True)
    parser.add_argument("--output-genome", type=Path, required=True)
    parser.add_argument("--output-annotation", type=Path, required=True)
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

    content = f"""# Organelle Submission Files

This directory contains files staged for submission of the {label} genome for `{args.assembly_name}` (version `{args.assembly_version}`).

## File Provenance

| Submission file | Source file | How the submission file was generated |
| --- | --- | --- |
| {quote_path(args.output_genome)} | {quote_path(args.input_genome)} | The organelle genome FASTA was copied without sequence changes and compressed with `gzip`. |
| {quote_path(args.output_annotation)} | {quote_path(args.input_annotation)} | The organelle annotation file emitted by Oatk was copied without content changes and compressed with `gzip`. |

## Notes

- These files were staged automatically from the Oatk outputs selected by `oatk_organelle` in `config/config.yml`.
- No automatic identifier renaming or annotation rewriting was applied to the organelle inputs.
- This README was generated automatically by `workflow/scripts/write_organelle_submission_readme.py`.
"""

    args.output_readme.write_text(content)


def main():
    args = parse_args()
    write_readme(args)


if __name__ == "__main__":
    main()

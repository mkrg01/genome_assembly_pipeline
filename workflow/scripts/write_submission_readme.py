import argparse
from pathlib import Path

from rename_submission_gene_models import make_submission_prefix


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write a README.md describing submission file provenance."
    )
    parser.add_argument("--assembly-name", required=True)
    parser.add_argument("--assembly-version", required=True)
    parser.add_argument("--input-assembly", type=Path, required=True)
    parser.add_argument("--input-isoform-cds", type=Path, required=True)
    parser.add_argument("--input-isoform-gff3", type=Path, required=True)
    parser.add_argument("--input-longest-cds", type=Path, required=True)
    parser.add_argument("--input-longest-gff3", type=Path, required=True)
    parser.add_argument("--output-assembly", type=Path, required=True)
    parser.add_argument("--output-isoform-cds", type=Path, required=True)
    parser.add_argument("--output-isoform-gff3", type=Path, required=True)
    parser.add_argument("--output-longest-cds", type=Path, required=True)
    parser.add_argument("--output-longest-gff3", type=Path, required=True)
    parser.add_argument("--output-readme", type=Path, required=True)
    return parser.parse_args()


def quote_path(path: Path):
    return f"`{path.as_posix()}`"


def write_readme(args):
    prefix = make_submission_prefix(args.assembly_name)
    args.output_readme.parent.mkdir(parents=True, exist_ok=True)

    content = f"""# Submission Files

This directory contains files prepared by the `format_for_submission` rule for the assembly `{args.assembly_name}` (version `{args.assembly_version}`).

## ID Naming Scheme

The submission annotation files use gene and transcript IDs derived from the assembly name `{args.assembly_name}`.

- Genus and species were inferred by replacing underscores with spaces.
- The prefix `{prefix}` was built from the first three letters of the genus and the first two letters of the species.
- Gene IDs were renamed from the BRAKER-style format `gN` to `{prefix}_000001`, `{prefix}_000002`, and so on.
- Transcript IDs retained their isoform suffixes. For example, `g1.t3` became `{prefix}_000001.t3`.
- In the GFF3 files, `ID` and `Parent` attributes were updated so that all parent-child relationships remained consistent after renaming.

## File Provenance

| Submission file | Source file | How the submission file was generated |
| --- | --- | --- |
| {quote_path(args.output_assembly)} | {quote_path(args.input_assembly)} | The RepeatMasker-masked genome FASTA was copied without sequence changes and compressed with `gzip`. |
| {quote_path(args.output_isoform_cds)} | {quote_path(args.input_isoform_cds)} | The CDS sequences were kept unchanged, but each FASTA header was renamed from the original BRAKER transcript ID to the submission ID scheme described above. The result was then compressed with `gzip`. |
| {quote_path(args.output_isoform_gff3)} | {quote_path(args.input_isoform_gff3)} | The original BRAKER GFF3 was rewritten to match the submission ID scheme. Gene `ID` values, mRNA `ID` values, child feature `ID` values, and all `Parent` attributes were updated accordingly. The result was then compressed with `gzip`. |
| {quote_path(args.output_longest_cds)} | {quote_path(args.input_longest_cds)} | The representative CDS FASTA, which had already been reduced upstream to one representative transcript per gene, was renamed to the submission ID scheme at the FASTA-header level and compressed with `gzip`. |
| {quote_path(args.output_longest_gff3)} | {quote_path(args.input_longest_gff3)} | The representative-model GFF3, which had already been reduced upstream to the representative transcript set, was rewritten so that all `ID` and `Parent` values matched the submission ID scheme, then compressed with `gzip`. |

## Notes

- The submission-ready annotation files were generated with `workflow/scripts/rename_submission_gene_models.py`.
- This README was generated automatically by `workflow/scripts/write_submission_readme.py`.
"""

    args.output_readme.write_text(content)


def main():
    args = parse_args()
    write_readme(args)


if __name__ == "__main__":
    main()

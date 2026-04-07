import argparse
import gzip
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter a GFF3 file to retain only genes, mRNAs, and mRNA child features "
            "associated with transcript IDs present in a FASTA file. Output order is "
            "preserved from the original GFF3."
        )
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="FASTA file whose record IDs define the mRNAs to retain.",
    )
    parser.add_argument(
        "--gff3",
        type=Path,
        required=True,
        help="Input GFF3 file to filter.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output GFF3 file.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_fasta_ids(path: Path):
    ids = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue

            record_id = line[1:].strip().split(None, 1)[0]
            if not record_id:
                raise ValueError(f"Missing FASTA record ID at line {line_number} in {path}")
            ids.add(record_id)

    if not ids:
        raise ValueError(f"No FASTA record IDs were found in {path}")

    return ids


def parse_attributes(attribute_text: str):
    attributes = {}
    for field in attribute_text.strip().split(";"):
        if not field or "=" not in field:
            continue
        key, value = field.split("=", 1)
        attributes[key] = value
    return attributes


def split_parents(parent_value: str):
    if not parent_value:
        return []
    return [parent for parent in parent_value.split(",") if parent]


def iter_gff_records(path: Path):
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("#") or not raw_line.strip():
                yield line_number, raw_line, None
                continue

            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Expected 9 columns in GFF3 at line {line_number} of {path}, "
                    f"but found {len(fields)}"
                )

            yield line_number, raw_line, fields


def collect_gene_ids(gff3_path: Path, keep_mrna_ids):
    keep_gene_ids = set()
    found_mrna_ids = set()

    for _, _, fields in iter_gff_records(gff3_path):
        if fields is None or fields[2] != "mRNA":
            continue

        attributes = parse_attributes(fields[8])
        mrna_id = attributes.get("ID")
        if mrna_id not in keep_mrna_ids:
            continue

        found_mrna_ids.add(mrna_id)
        keep_gene_ids.update(split_parents(attributes.get("Parent", "")))

    return keep_gene_ids, found_mrna_ids


def write_filtered_gff3(gff3_path: Path, output_path: Path, keep_mrna_ids, keep_gene_ids):
    kept_feature_lines = 0
    kept_gene_lines = 0
    kept_mrna_lines = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as out_handle:
        for _, raw_line, fields in iter_gff_records(gff3_path):
            if fields is None:
                out_handle.write(raw_line)
                continue

            feature_type = fields[2]
            attributes = parse_attributes(fields[8])
            feature_id = attributes.get("ID", "")
            parent_ids = set(split_parents(attributes.get("Parent", "")))

            if feature_type == "gene" and feature_id in keep_gene_ids:
                kept_gene_lines += 1
                kept_feature_lines += 1
                out_handle.write(raw_line)
                continue

            if feature_type == "mRNA" and feature_id in keep_mrna_ids:
                kept_mrna_lines += 1
                kept_feature_lines += 1
                out_handle.write(raw_line)
                continue

            if parent_ids & keep_mrna_ids:
                kept_feature_lines += 1
                out_handle.write(raw_line)

    return kept_gene_lines, kept_mrna_lines, kept_feature_lines


def main():
    args = parse_args()
    keep_mrna_ids = parse_fasta_ids(args.fasta)
    keep_gene_ids, found_mrna_ids = collect_gene_ids(args.gff3, keep_mrna_ids)

    missing_mrna_ids = sorted(keep_mrna_ids - found_mrna_ids)
    if missing_mrna_ids:
        preview = ", ".join(missing_mrna_ids[:10])
        suffix = "" if len(missing_mrna_ids) <= 10 else ", ..."
        print(
            f"Warning: {len(missing_mrna_ids)} FASTA IDs were not found as mRNA IDs in "
            f"{args.gff3}: {preview}{suffix}",
            file=sys.stderr,
        )

    kept_gene_lines, kept_mrna_lines, kept_feature_lines = write_filtered_gff3(
        args.gff3, args.output, found_mrna_ids, keep_gene_ids
    )

    print(
        f"Wrote {args.output} with {kept_gene_lines} gene lines, "
        f"{kept_mrna_lines} mRNA lines, and {kept_feature_lines} total feature lines.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

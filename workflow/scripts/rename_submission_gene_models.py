import argparse
import gzip
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rename BRAKER gene and transcript IDs for submission outputs using a "
            "prefix derived from the species name, and apply the same mapping to "
            "FASTA and GFF3 files."
        )
    )
    parser.add_argument(
        "--species-name",
        required=True,
        help='Scientific name used to derive the prefix (e.g., "Plumbago_auriculata").',
    )
    parser.add_argument("--isoform-cds", type=Path, required=True)
    parser.add_argument("--isoform-gff3", type=Path, required=True)
    parser.add_argument("--longest-cds", type=Path, required=True)
    parser.add_argument("--longest-gff3", type=Path, required=True)
    parser.add_argument("--output-isoform-cds", type=Path, required=True)
    parser.add_argument("--output-isoform-gff3", type=Path, required=True)
    parser.add_argument("--output-longest-cds", type=Path, required=True)
    parser.add_argument("--output-longest-gff3", type=Path, required=True)
    return parser.parse_args()


def open_text_read(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def open_text_write(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt")
    return path.open("w")


def split_parents(parent_value: str):
    if not parent_value:
        return []
    return [parent for parent in parent_value.split(",") if parent]


def parse_attributes(attribute_text: str):
    attributes = {}
    for field in attribute_text.strip().split(";"):
        if not field or "=" not in field:
            continue
        key, value = field.split("=", 1)
        attributes[key] = value
    return attributes


def format_attributes(attributes):
    if not attributes:
        return "."
    return ";".join(f"{key}={value}" for key, value in attributes.items()) + ";"


def make_submission_prefix(species_name: str):
    tokens = species_name.replace("_", " ").split()
    if len(tokens) < 2:
        raise ValueError(
            "Species name must contain at least genus and species, "
            f"but got: {species_name!r}"
        )

    genus = tokens[0]
    species = tokens[1]
    return f"{genus[:3].capitalize()}{species[:2].lower()}"


def iter_gff_records(path: Path):
    with open_text_read(path) as handle:
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


def build_transcript_id(old_gene_id: str, new_gene_id: str, transcript_id: str):
    if not transcript_id.startswith(old_gene_id):
        raise ValueError(
            f"Transcript ID {transcript_id!r} does not start with its parent gene ID "
            f"{old_gene_id!r}"
        )
    return f"{new_gene_id}{transcript_id[len(old_gene_id):]}"


def build_identifier_maps(gff3_path: Path, prefix: str):
    gene_map = {}
    transcript_map = {}
    gene_counter = 0

    for line_number, _, fields in iter_gff_records(gff3_path):
        if fields is None:
            continue

        feature_type = fields[2]
        attributes = parse_attributes(fields[8])

        if feature_type == "gene":
            gene_id = attributes.get("ID")
            if not gene_id:
                raise ValueError(
                    f"Missing gene ID at line {line_number} of {gff3_path}"
                )
            if gene_id not in gene_map:
                gene_counter += 1
                gene_map[gene_id] = f"{prefix}_{gene_counter:06d}"
            continue

        if feature_type != "mRNA":
            continue

        transcript_id = attributes.get("ID")
        parent_ids = split_parents(attributes.get("Parent", ""))
        if not transcript_id:
            raise ValueError(
                f"Missing mRNA ID at line {line_number} of {gff3_path}"
            )
        if len(parent_ids) != 1:
            raise ValueError(
                f"Expected exactly one parent gene for mRNA {transcript_id!r} at "
                f"line {line_number} of {gff3_path}, but found {parent_ids}"
            )

        gene_id = parent_ids[0]
        if gene_id not in gene_map:
            raise ValueError(
                f"Parent gene {gene_id!r} for transcript {transcript_id!r} was not "
                f"seen before line {line_number} of {gff3_path}"
            )

        transcript_map[transcript_id] = build_transcript_id(
            gene_id, gene_map[gene_id], transcript_id
        )

    if not gene_map:
        raise ValueError(f"No gene features were found in {gff3_path}")
    if not transcript_map:
        raise ValueError(f"No mRNA features were found in {gff3_path}")

    return gene_map, transcript_map


def rename_parent_id(parent_id: str, gene_map, transcript_map):
    if parent_id in transcript_map:
        return transcript_map[parent_id]
    if parent_id in gene_map:
        return gene_map[parent_id]
    raise ValueError(f"Could not rename parent ID {parent_id!r}")


def rename_child_feature_id(feature_id: str, old_parents, new_parents, gene_map, transcript_map):
    if feature_id in transcript_map:
        return transcript_map[feature_id]
    if feature_id in gene_map:
        return gene_map[feature_id]

    for old_parent, new_parent in zip(old_parents, new_parents):
        if feature_id.startswith(f"{old_parent}."):
            return f"{new_parent}{feature_id[len(old_parent):]}"

    return feature_id


def rewrite_gff3(input_path: Path, output_path: Path, gene_map, transcript_map):
    feature_count = 0

    with open_text_write(output_path) as out_handle:
        for line_number, raw_line, fields in iter_gff_records(input_path):
            if fields is None:
                out_handle.write(raw_line)
                continue

            feature_type = fields[2]
            attributes = parse_attributes(fields[8])

            if feature_type == "gene":
                gene_id = attributes.get("ID")
                if gene_id not in gene_map:
                    raise ValueError(
                        f"Unknown gene ID {gene_id!r} at line {line_number} of {input_path}"
                    )
                attributes["ID"] = gene_map[gene_id]
            elif feature_type == "mRNA":
                transcript_id = attributes.get("ID")
                parent_ids = split_parents(attributes.get("Parent", ""))
                if transcript_id not in transcript_map:
                    raise ValueError(
                        f"Unknown transcript ID {transcript_id!r} at line {line_number} "
                        f"of {input_path}"
                    )
                attributes["ID"] = transcript_map[transcript_id]
                attributes["Parent"] = ",".join(
                    rename_parent_id(parent_id, gene_map, transcript_map)
                    for parent_id in parent_ids
                )
            else:
                old_parents = split_parents(attributes.get("Parent", ""))
                if old_parents:
                    new_parents = [
                        rename_parent_id(parent_id, gene_map, transcript_map)
                        for parent_id in old_parents
                    ]
                    attributes["Parent"] = ",".join(new_parents)
                    if "ID" in attributes:
                        attributes["ID"] = rename_child_feature_id(
                            attributes["ID"],
                            old_parents,
                            new_parents,
                            gene_map,
                            transcript_map,
                        )

            fields[8] = format_attributes(attributes)
            out_handle.write("\t".join(fields) + "\n")
            feature_count += 1

    return feature_count


def rewrite_fasta(input_path: Path, output_path: Path, transcript_map):
    record_count = 0

    with open_text_read(input_path) as in_handle, open_text_write(output_path) as out_handle:
        for line_number, line in enumerate(in_handle, start=1):
            if not line.startswith(">"):
                out_handle.write(line)
                continue

            header = line[1:].rstrip("\n")
            if not header:
                raise ValueError(f"Missing FASTA header at line {line_number} of {input_path}")

            header_fields = header.split(None, 1)
            transcript_id = header_fields[0]
            description = f" {header_fields[1]}" if len(header_fields) == 2 else ""

            if transcript_id not in transcript_map:
                raise ValueError(
                    f"FASTA ID {transcript_id!r} at line {line_number} of {input_path} "
                    "was not found in the transcript mapping"
                )

            out_handle.write(f">{transcript_map[transcript_id]}{description}\n")
            record_count += 1

    if record_count == 0:
        raise ValueError(f"No FASTA records were found in {input_path}")

    return record_count


def main():
    args = parse_args()
    prefix = make_submission_prefix(args.species_name)
    gene_map, transcript_map = build_identifier_maps(args.isoform_gff3, prefix)

    isoform_record_count = rewrite_fasta(
        args.isoform_cds, args.output_isoform_cds, transcript_map
    )
    longest_record_count = rewrite_fasta(
        args.longest_cds, args.output_longest_cds, transcript_map
    )
    isoform_feature_count = rewrite_gff3(
        args.isoform_gff3, args.output_isoform_gff3, gene_map, transcript_map
    )
    longest_feature_count = rewrite_gff3(
        args.longest_gff3, args.output_longest_gff3, gene_map, transcript_map
    )

    print(
        f"Submission gene prefix: {prefix}",
        file=sys.stderr,
    )
    print(
        f"Renamed {len(gene_map)} genes and {len(transcript_map)} transcripts.",
        file=sys.stderr,
    )
    print(
        f"Wrote {args.output_isoform_cds} ({isoform_record_count} records) and "
        f"{args.output_longest_cds} ({longest_record_count} records).",
        file=sys.stderr,
    )
    print(
        f"Wrote {args.output_isoform_gff3} ({isoform_feature_count} features) and "
        f"{args.output_longest_gff3} ({longest_feature_count} features).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

import json
import re
import shutil
from pathlib import Path


GENBANK_SUFFIXES = (".gb", ".gbk", ".gbf")
FEATURE_INDENT = "     "
QUALIFIER_INDENT = "                     "
SPECIES_PLACEHOLDER = "YOUR_SPECIES"
GENBANK_TOPOLOGIES = {"circular", "linear"}
FASTA_CIRCULAR_PATTERN = re.compile(r"(?:^|\s)circular=(true|false)(?=\s|$)", re.I)
GENBANK_DATE_PATTERN = re.compile(r"\d{2}-[A-Z]{3}-\d{4}")
GENBANK_MOLECULE_PATTERN = re.compile(r"\b(mRNA|DNA|RNA)\b")


def organism_name_from_assembly_name(assembly_name: str):
    organism_name = " ".join(assembly_name.replace("_", " ").split())
    if not organism_name:
        raise ValueError("Assembly name must not be empty when deriving organism name.")
    if any(character in organism_name for character in '"\r\n'):
        raise ValueError(
            f"Assembly name cannot be used as a GenBank organism: {assembly_name!r}"
        )
    return organism_name


def copy_first_genbank(output_dir: Path, annotation_path: Path):
    candidates = sorted(
        path
        for path in output_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in GENBANK_SUFFIXES
    )
    if not candidates:
        raise RuntimeError(f"No GenBank output was found under {output_dir}.")

    annotation_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(candidates[0], annotation_path)
    return candidates[0]


def split_genbank_records(lines):
    records = []
    current = []
    for line in lines:
        current.append(line)
        if line.strip() == "//":
            records.append(current)
            current = []
    if current:
        records.append(current)
    return records


def trim_genbank_record_to_core_sections(record, path: Path):
    locus_index = None
    features_index = None
    origin_index = None
    for index, line in enumerate(record):
        if line.startswith("LOCUS") and locus_index is None:
            locus_index = index
        elif line.startswith("FEATURES") and features_index is None:
            features_index = index
        elif line.startswith("ORIGIN") and origin_index is None:
            origin_index = index

    missing = [
        name
        for name, index in (
            ("LOCUS", locus_index),
            ("FEATURES", features_index),
            ("ORIGIN", origin_index),
        )
        if index is None
    ]
    if missing:
        raise ValueError(
            f"GenBank record in {path} is missing required section(s): "
            + ", ".join(missing)
        )
    if not locus_index < features_index < origin_index:
        raise ValueError(
            f"GenBank record in {path} has sections out of order: "
            "expected LOCUS, FEATURES, ORIGIN."
        )

    return [record[locus_index], *record[features_index:]]


def trim_genbank_to_core_sections(path: Path):
    lines = path.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    trimmed_records = [
        trim_genbank_record_to_core_sections(record, path) for record in records
    ]
    trimmed_lines = [line for record in trimmed_records for line in record]
    changed = lines != trimmed_lines
    if changed:
        path.write_text("".join(trimmed_lines))
    return {
        "core_sections": ["LOCUS", "FEATURES", "ORIGIN"],
        "core_sections_record_count": len(trimmed_records),
        "core_sections_changed": changed,
    }


def topology_from_fasta_header(header: str):
    match = FASTA_CIRCULAR_PATTERN.search(header)
    if match is None:
        return None
    return "circular" if match.group(1).lower() == "true" else "linear"


def normalize_genbank_topology(topology):
    if topology is None:
        return None
    normalized = str(topology).strip().lower()
    if normalized not in GENBANK_TOPOLOGIES:
        raise ValueError(
            f"Unsupported GenBank topology '{topology}'. Expected 'circular' or 'linear'."
        )
    return normalized


def extract_locus_metadata(line: str, sequence_length=None):
    if not line.startswith("LOCUS"):
        raise ValueError(f"Expected a GenBank LOCUS line, got: {line!r}")

    length_match = re.search(r"(\d+)\s+bp\b", line)
    if length_match:
        parsed_length = int(length_match.group(1))
        name = line[len("LOCUS") : length_match.start(1)].strip()
        suffix = line[length_match.end() :]
    else:
        parsed_length = None
        name = line[len("LOCUS") :].strip()
        suffix = ""

    length = sequence_length if sequence_length is not None else parsed_length
    if sequence_length is not None:
        length_text = str(sequence_length)
        if name.endswith(length_text):
            name = name[: -len(length_text)].rstrip()

    molecule_match = GENBANK_MOLECULE_PATTERN.search(suffix)
    molecule = molecule_match.group(1) if molecule_match else "DNA"

    topology_match = re.search(r"circular|linear", suffix, flags=re.I)
    topology = topology_match.group(0).lower() if topology_match else None

    date_match = GENBANK_DATE_PATTERN.search(suffix)
    date = date_match.group(0) if date_match else None

    division = None
    if date_match:
        before_date = suffix[: date_match.start()]
        division_match = re.search(
            r"(?:circular|linear)?\s*([A-Z]{3})\s*$",
            before_date,
            flags=re.I,
        )
        if division_match:
            division = division_match.group(1).upper()
    if division is None:
        division = "PLN"

    return {
        "name": name,
        "length": length,
        "molecule": molecule,
        "topology": topology,
        "division": division,
        "date": date,
    }


def format_locus_line(
    original_line: str,
    *,
    topology=None,
    locus_name=None,
    sequence_length=None,
):
    metadata = extract_locus_metadata(original_line, sequence_length=sequence_length)
    topology = normalize_genbank_topology(topology) or metadata["topology"]
    name = locus_name or metadata["name"]
    length = sequence_length if sequence_length is not None else metadata["length"]
    if length is None:
        raise ValueError(f"Could not determine LOCUS sequence length from: {original_line!r}")

    if topology is None:
        topology_text = ""
    else:
        topology_text = f"{topology:<8} "
    trailing = " ".join(item for item in (metadata["division"], metadata["date"]) if item)
    return (
        f"LOCUS       {name:<24} {length:>11} bp    "
        f"{metadata['molecule']:<6} {topology_text}{trailing}"
    ).rstrip()


def replace_locus_topology(line: str, topology: str):
    return format_locus_line(line, topology=topology)


def curate_genbank_locus(
    path: Path,
    *,
    topology=None,
    locus_name=None,
    sequence_length=None,
):
    topology = normalize_genbank_topology(topology)
    if topology is None and locus_name is None and sequence_length is None:
        return {
            "locus_topology_before": None,
            "locus_topology_after": None,
            "locus_topology_changed": False,
            "locus_topology_skipped_reason": "no input FASTA circular flag provided",
            "locus_line_before": None,
            "locus_line_after": None,
            "locus_line_changed": False,
        }

    lines = path.read_text().splitlines(keepends=True)
    locus_index = None
    for index, line in enumerate(lines):
        if line.startswith("LOCUS"):
            locus_index = index
            break
    if locus_index is None:
        return {
            "locus_topology_before": None,
            "locus_topology_after": topology,
            "locus_topology_changed": False,
            "locus_topology_skipped_reason": "no LOCUS line found",
            "locus_line_before": None,
            "locus_line_after": None,
            "locus_line_changed": False,
        }

    line = lines[locus_index]
    newline = "\n" if line.endswith("\n") else ""
    line_without_newline = line.rstrip("\n")
    before_metadata = extract_locus_metadata(
        line_without_newline,
        sequence_length=sequence_length,
    )
    before = before_metadata["topology"]
    updated_line_without_newline = format_locus_line(
        line_without_newline,
        topology=topology,
        locus_name=locus_name,
        sequence_length=sequence_length,
    )
    updated_line = updated_line_without_newline + newline
    changed = updated_line != line
    if changed:
        lines[locus_index] = updated_line
        path.write_text("".join(lines))

    return {
        "locus_topology_before": before,
        "locus_topology_after": topology or before,
        "locus_topology_changed": before != (topology or before),
        "locus_topology_skipped_reason": None,
        "locus_line_before": line_without_newline,
        "locus_line_after": updated_line_without_newline,
        "locus_line_changed": changed,
    }


def curate_genbank_locus_topology(path: Path, topology):
    return curate_genbank_locus(path, topology=topology)


def curate_source_organism(path: Path, assembly_name: str):
    organism_name = organism_name_from_assembly_name(assembly_name)
    organism_line = f'{QUALIFIER_INDENT}/organism="{organism_name}"\n'
    lines = path.read_text().splitlines(keepends=True)

    source_index = None
    for index, line in enumerate(lines):
        if line.startswith(f"{FEATURE_INDENT}source"):
            source_index = index
            break
    if source_index is None:
        raise ValueError(f"No source feature was found in {path}.")

    next_feature_index = len(lines)
    for index in range(source_index + 1, len(lines)):
        line = lines[index]
        if line.startswith(FEATURE_INDENT) and not line.startswith(QUALIFIER_INDENT):
            next_feature_index = index
            break

    previous_organism = None
    organism_index = None
    for index in range(source_index + 1, next_feature_index):
        stripped = lines[index].strip()
        if stripped.startswith("/organism="):
            organism_index = index
            previous_organism = stripped.removeprefix("/organism=").strip('"')
            break

    if organism_index is None:
        lines.insert(source_index + 1, organism_line)
    elif lines[organism_index] != organism_line:
        lines[organism_index] = organism_line

    changed = previous_organism != organism_name
    if changed:
        path.write_text("".join(lines))

    return {
        "source_organism_before": previous_organism,
        "source_organism_after": organism_name,
        "source_organism_changed": changed,
    }


def curate_source_organelle(path: Path, organelle: str | None = None):
    if organelle is None:
        return {
            "source_organelle_before": None,
            "source_organelle_after": None,
            "source_organelle_changed": False,
            "source_organelle_skipped_reason": "no organelle qualifier requested",
        }
    if any(character in organelle for character in '"\r\n'):
        raise ValueError(f"Invalid source organelle qualifier: {organelle!r}")

    organelle_line = f'{QUALIFIER_INDENT}/organelle="{organelle}"\n'
    lines = path.read_text().splitlines(keepends=True)

    source_index = None
    for index, line in enumerate(lines):
        if line.startswith(f"{FEATURE_INDENT}source"):
            source_index = index
            break
    if source_index is None:
        raise ValueError(f"No source feature was found in {path}.")

    next_feature_index = len(lines)
    for index in range(source_index + 1, len(lines)):
        line = lines[index]
        if line.startswith(FEATURE_INDENT) and not line.startswith(QUALIFIER_INDENT):
            next_feature_index = index
            break

    previous_organelle = None
    organelle_index = None
    organism_index = None
    for index in range(source_index + 1, next_feature_index):
        stripped = lines[index].strip()
        if stripped.startswith("/organism="):
            organism_index = index
        elif stripped.startswith("/organelle="):
            organelle_index = index
            previous_organelle = stripped.removeprefix("/organelle=").strip('"')
            break

    if organelle_index is None:
        insert_index = organism_index + 1 if organism_index is not None else source_index + 1
        lines.insert(insert_index, organelle_line)
    elif lines[organelle_index] != organelle_line:
        lines[organelle_index] = organelle_line

    changed = previous_organelle != organelle
    if changed:
        path.write_text("".join(lines))

    return {
        "source_organelle_before": previous_organelle,
        "source_organelle_after": organelle,
        "source_organelle_changed": changed,
        "source_organelle_skipped_reason": None,
    }


def normalize_taxid(taxid):
    if taxid is None:
        return None
    if isinstance(taxid, bool):
        raise ValueError("taxid must not be a boolean.")
    value = str(taxid).strip()
    if not value:
        return None
    if not value.isdigit() or int(value) < 1:
        raise ValueError("taxid must be a positive integer NCBI Taxonomy ID.")
    return value


def curate_source_taxon_db_xref(path: Path, taxid=None, remove_without_taxid=False):
    taxid = normalize_taxid(taxid)
    desired_xref = f"taxon:{taxid}" if taxid else None
    previous_xrefs = []
    lines = path.read_text().splitlines(keepends=True)

    source_index = None
    for index, line in enumerate(lines):
        if line.startswith(f"{FEATURE_INDENT}source"):
            source_index = index
            break
    if source_index is None:
        raise ValueError(f"No source feature was found in {path}.")

    next_feature_index = len(lines)
    for index in range(source_index + 1, len(lines)):
        line = lines[index]
        if line.startswith(FEATURE_INDENT) and not line.startswith(QUALIFIER_INDENT):
            next_feature_index = index
            break

    curated_lines = []
    organism_curated_index = None
    first_taxon_curated_index = None
    for index, line in enumerate(lines):
        stripped = line.strip()
        in_source = source_index < index < next_feature_index
        if in_source and stripped.startswith('/db_xref="taxon:'):
            if first_taxon_curated_index is None:
                first_taxon_curated_index = len(curated_lines)
            previous_xrefs.append(stripped.removeprefix("/db_xref=").strip('"'))
            continue
        if in_source and stripped.startswith("/organism="):
            organism_curated_index = len(curated_lines)
        curated_lines.append(line)

    if desired_xref is not None:
        xref_line = f'{QUALIFIER_INDENT}/db_xref="{desired_xref}"\n'
        if first_taxon_curated_index is not None:
            insert_index = first_taxon_curated_index
        elif organism_curated_index is not None:
            insert_index = organism_curated_index + 1
        else:
            insert_index = source_index + 1
        curated_lines.insert(insert_index, xref_line)
        current_xrefs = [desired_xref]
    elif previous_xrefs and not remove_without_taxid:
        insert_index = first_taxon_curated_index
        for xref in previous_xrefs:
            curated_lines.insert(
                insert_index,
                f'{QUALIFIER_INDENT}/db_xref="{xref}"\n',
            )
            insert_index += 1
        current_xrefs = previous_xrefs
    else:
        current_xrefs = []

    changed = previous_xrefs != current_xrefs
    if changed:
        path.write_text("".join(curated_lines))

    return {
        "source_taxon_db_xrefs_before": previous_xrefs,
        "source_taxon_db_xrefs_after": current_xrefs,
        "source_taxon_db_xref_changed": changed,
    }


def curate_species_placeholders(path: Path, assembly_name: str):
    organism_name = organism_name_from_assembly_name(assembly_name)
    lines = path.read_text().splitlines(keepends=True)
    replacements = 0
    curated_lines = []
    for line in lines:
        replacements += line.count(SPECIES_PLACEHOLDER)
        curated_lines.append(line.replace(SPECIES_PLACEHOLDER, organism_name))

    changed = replacements > 0
    if changed:
        path.write_text("".join(curated_lines))

    return {
        "species_placeholder": SPECIES_PLACEHOLDER,
        "species_placeholder_replacement": organism_name,
        "species_placeholder_replacements": replacements,
        "species_placeholder_changed": changed,
    }


def curate_genbank_source_metadata(
    path: Path,
    assembly_name: str,
    taxid=None,
    organelle=None,
):
    normalized_taxid = normalize_taxid(taxid)
    source_organism = curate_source_organism(path, assembly_name)
    source_organelle = curate_source_organelle(path, organelle)
    source_taxon_db_xref = curate_source_taxon_db_xref(
        path,
        taxid=normalized_taxid,
        remove_without_taxid=source_organism["source_organism_changed"],
    )
    species_placeholders = curate_species_placeholders(path, assembly_name)
    return {
        "source_organism": source_organism,
        "source_organelle": source_organelle,
        "source_taxon_db_xref": source_taxon_db_xref,
        "species_placeholders": species_placeholders,
    }


def append_post_curation_summary(lines, post_curation):
    core_sections = post_curation.get("core_sections")
    if core_sections:
        sections = ", ".join(core_sections["core_sections"])
        record_count = core_sections["core_sections_record_count"]
        if core_sections["core_sections_changed"]:
            lines.append(
                f"- GenBank file: trimmed {record_count} record(s) to {sections}"
            )
        else:
            lines.append(
                f"- GenBank file: already contained only {sections} section(s)"
            )

    locus_topology = post_curation.get("locus_topology")
    if locus_topology:
        before = locus_topology["locus_topology_before"]
        after = locus_topology["locus_topology_after"]
        skipped_reason = locus_topology["locus_topology_skipped_reason"]
        if skipped_reason:
            lines.append(f"- LOCUS topology: skipped ({skipped_reason})")
        elif before is None:
            lines.append(f'- LOCUS topology: added "{after}"')
        elif locus_topology["locus_topology_changed"]:
            lines.append(f'- LOCUS topology: changed from "{before}" to "{after}"')
        else:
            lines.append(f'- LOCUS topology: already "{after}"')

    source_organism = post_curation.get("source_organism")
    if source_organism:
        before = source_organism["source_organism_before"]
        after = source_organism["source_organism_after"]
        if before is None:
            lines.append(f'- source /organism: added "{after}"')
        elif source_organism["source_organism_changed"]:
            lines.append(f'- source /organism: changed from "{before}" to "{after}"')
        else:
            lines.append(f'- source /organism: already "{after}"')

    source_organelle = post_curation.get("source_organelle")
    if source_organelle and not source_organelle["source_organelle_skipped_reason"]:
        before = source_organelle["source_organelle_before"]
        after = source_organelle["source_organelle_after"]
        if before is None:
            lines.append(f'- source /organelle: added "{after}"')
        elif source_organelle["source_organelle_changed"]:
            lines.append(f'- source /organelle: changed from "{before}" to "{after}"')
        else:
            lines.append(f'- source /organelle: already "{after}"')

    source_taxon_db_xref = post_curation.get("source_taxon_db_xref")
    if source_taxon_db_xref and source_taxon_db_xref["source_taxon_db_xref_changed"]:
        before = source_taxon_db_xref["source_taxon_db_xrefs_before"]
        after = source_taxon_db_xref["source_taxon_db_xrefs_after"]
        before_text = ", ".join(before) if before else "none"
        after_text = ", ".join(after) if after else "none"
        lines.append(
            f"- source taxon db_xref: changed from {before_text} to {after_text}"
        )

    species_placeholders = post_curation.get("species_placeholders")
    if species_placeholders and species_placeholders["species_placeholder_changed"]:
        placeholder = species_placeholders["species_placeholder"]
        replacement = species_placeholders["species_placeholder_replacement"]
        count = species_placeholders["species_placeholder_replacements"]
        lines.append(
            f'- species placeholders: replaced {count} "{placeholder}" '
            f'occurrence(s) with "{replacement}"'
        )


def format_post_curation_record(post_curation, annotation_path: Path | None = None):
    if annotation_path is None:
        manual_note_line = (
            "This file is regenerated by the workflow. If you edit the GenBank "
            "annotation after the organelle annotation step finishes, record those "
            "edits below."
        )
    else:
        manual_note_line = (
            "This file is regenerated by the workflow. If you edit "
            f"`{annotation_path.name}` after the organelle annotation step "
            "finishes, record those edits below."
        )

    lines = [
        manual_note_line,
        "",
        "# Automatic post-curation",
        "",
    ]

    record_curations = post_curation.get("records")
    if record_curations:
        lines.append(
            "- PMGA was run separately for each input FASTA record; the curated "
            "GenBank records were then concatenated into the final multi-record file."
        )
        for record_curation in record_curations:
            record_id = record_curation.get("input_record_id")
            record_annotation = record_curation.get("annotation")
            label = f"`{record_id}`" if record_id else "record"
            if record_annotation:
                label = f"{label} -> `{Path(record_annotation).name}`"
            lines.append(f"- {label}")
        lines.append("")

        for record_curation in record_curations:
            record_id = record_curation.get("input_record_id", "record")
            lines.extend([f"## {record_id}", ""])
            append_post_curation_summary(
                lines,
                record_curation.get("post_curation", {}),
            )
            lines.append("")
    else:
        append_post_curation_summary(lines, post_curation)

    lines.extend(
        [
            "",
            "# Manual post-curation",
            "",
        ]
    )
    return "\n".join(lines)


def write_post_curation_record(
    path: Path,
    post_curation,
    annotation_path: Path | None = None,
):
    updated = format_post_curation_record(post_curation, annotation_path)

    path.parent.mkdir(parents=True, exist_ok=True)
    changed = (not path.exists()) or path.read_text() != updated
    if changed:
        path.write_text(updated)

    return {
        "post_curation_record": str(path),
        "post_curation_record_changed": changed,
    }


def write_run_manifest(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")

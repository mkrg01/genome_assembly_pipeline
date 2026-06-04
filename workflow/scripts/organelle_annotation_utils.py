import json
import re
import shutil
from pathlib import Path


GENBANK_SUFFIXES = (".gb", ".gbk", ".gbf")
FEATURE_INDENT = "     "
QUALIFIER_INDENT = "                     "
ORGANISM_PREFIX = "  ORGANISM  "
ORGANISM_TAXONOMY_INDENT = "            "
GENBANK_LINE_WIDTH = 79
SPECIES_PLACEHOLDER = "YOUR_SPECIES"
OMITTED_GENBANK_TAXONOMY_NAMES = {"root", "cellular organisms"}
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


def load_taxonomy_lineage_record(path):
    if path is None:
        return None
    if not str(path).strip():
        return None
    path = Path(path)
    return json.loads(path.read_text())


def genbank_taxonomy_text(taxonomy_names):
    if not taxonomy_names:
        return ""
    return "; ".join(taxonomy_names) + "."


def wrap_genbank_taxonomy(taxonomy_names):
    items = [
        f"{name};" if index < len(taxonomy_names) - 1 else f"{name}."
        for index, name in enumerate(taxonomy_names)
    ]
    lines = []
    current = ""
    content_width = GENBANK_LINE_WIDTH - len(ORGANISM_TAXONOMY_INDENT)
    for item in items:
        candidate = item if not current else f"{current} {item}"
        if current and len(candidate) > content_width:
            lines.append(f"{ORGANISM_TAXONOMY_INDENT}{current}\n")
            current = item
        else:
            current = candidate
    if current:
        lines.append(f"{ORGANISM_TAXONOMY_INDENT}{current}\n")
    return lines


def taxonomy_record_from_ncbi_taxa(ncbi, taxid):
    taxid = normalize_taxid(taxid)
    lineage_taxids = [
        int(lineage_taxid) for lineage_taxid in ncbi.get_lineage(int(taxid))
    ]
    names_by_taxid = ncbi.get_taxid_translator(lineage_taxids)
    ranks_by_taxid = ncbi.get_rank(lineage_taxids)
    lineage = []
    genbank_taxonomy = []

    for lineage_taxid in lineage_taxids:
        name = names_by_taxid.get(lineage_taxid)
        if name is None:
            name = str(lineage_taxid)
        rank = ranks_by_taxid.get(lineage_taxid)
        lineage.append(
            {
                "taxid": str(lineage_taxid),
                "name": name,
                "rank": rank,
            }
        )
        if lineage_taxid == int(taxid):
            continue
        if name in OMITTED_GENBANK_TAXONOMY_NAMES:
            continue
        genbank_taxonomy.append(name)

    return {
        "source": "NCBI Taxonomy via ete4",
        "taxid": taxid,
        "scientific_name": names_by_taxid.get(int(taxid)),
        "lineage": lineage,
        "genbank_taxonomy": genbank_taxonomy,
        "genbank_taxonomy_text": genbank_taxonomy_text(genbank_taxonomy),
    }


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


def extract_taxonomy_names(taxonomy_lineage):
    if taxonomy_lineage is None:
        return []
    if isinstance(taxonomy_lineage, dict):
        taxonomy_names = taxonomy_lineage.get("genbank_taxonomy", [])
    else:
        taxonomy_names = taxonomy_lineage
    if taxonomy_names is None:
        return []
    normalized = []
    for name in taxonomy_names:
        if not isinstance(name, str):
            raise ValueError(f"Taxonomy lineage names must be strings: {name!r}")
        name = " ".join(name.split())
        if any(character in name for character in "\r\n;"):
            raise ValueError(f"Invalid taxonomy lineage name: {name!r}")
        if name:
            normalized.append(name)
    return normalized


def curate_organism_taxonomy_lineage(path: Path, taxonomy_lineage=None):
    taxonomy_names = extract_taxonomy_names(taxonomy_lineage)
    source = None
    taxid = None
    if isinstance(taxonomy_lineage, dict):
        source = taxonomy_lineage.get("source")
        taxid = taxonomy_lineage.get("taxid")

    if not taxonomy_names:
        return {
            "organism_taxonomy_source": source,
            "organism_taxonomy_taxid": taxid,
            "organism_taxonomy_before": None,
            "organism_taxonomy_after": None,
            "organism_taxonomy_changed": False,
            "organism_taxonomy_skipped_reason": "no taxonomy lineage provided",
        }

    lines = path.read_text().splitlines(keepends=True)
    organism_index = None
    for index, line in enumerate(lines):
        if line.startswith(ORGANISM_PREFIX):
            organism_index = index
            break

    if organism_index is None:
        return {
            "organism_taxonomy_source": source,
            "organism_taxonomy_taxid": taxid,
            "organism_taxonomy_before": None,
            "organism_taxonomy_after": genbank_taxonomy_text(taxonomy_names),
            "organism_taxonomy_changed": False,
            "organism_taxonomy_skipped_reason": "no ORGANISM line found",
        }

    taxonomy_start = organism_index + 1
    taxonomy_end = taxonomy_start
    while taxonomy_end < len(lines) and lines[taxonomy_end].startswith(
        ORGANISM_TAXONOMY_INDENT
    ):
        taxonomy_end += 1

    previous_lines = lines[taxonomy_start:taxonomy_end]
    next_lines = wrap_genbank_taxonomy(taxonomy_names)
    changed = [line.rstrip("\n") for line in previous_lines] != [
        line.rstrip("\n") for line in next_lines
    ]
    if changed:
        lines[taxonomy_start:taxonomy_end] = next_lines
        path.write_text("".join(lines))

    before_text = " ".join(line.strip() for line in previous_lines) or None
    after_text = genbank_taxonomy_text(taxonomy_names)
    return {
        "organism_taxonomy_source": source,
        "organism_taxonomy_taxid": taxid,
        "organism_taxonomy_before": before_text,
        "organism_taxonomy_after": after_text,
        "organism_taxonomy_changed": changed,
        "organism_taxonomy_skipped_reason": None,
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


def curate_genbank_species_metadata(
    path: Path,
    assembly_name: str,
    taxid=None,
    taxonomy_lineage=None,
):
    normalized_taxid = normalize_taxid(taxid)
    if isinstance(taxonomy_lineage, dict) and taxonomy_lineage.get("taxid"):
        taxonomy_taxid = normalize_taxid(taxonomy_lineage["taxid"])
        if normalized_taxid and taxonomy_taxid != normalized_taxid:
            raise ValueError(
                "Taxonomy lineage taxid does not match configured taxid: "
                f"{taxonomy_taxid} != {normalized_taxid}"
            )
    source_organism = curate_source_organism(path, assembly_name)
    source_taxon_db_xref = curate_source_taxon_db_xref(
        path,
        taxid=normalized_taxid,
        remove_without_taxid=source_organism["source_organism_changed"],
    )
    species_placeholders = curate_species_placeholders(path, assembly_name)
    organism_taxonomy = curate_organism_taxonomy_lineage(path, taxonomy_lineage)
    return {
        "source_organism": source_organism,
        "source_taxon_db_xref": source_taxon_db_xref,
        "species_placeholders": species_placeholders,
        "organism_taxonomy": organism_taxonomy,
    }


def append_post_curation_summary(lines, post_curation):
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

    source_taxon_db_xref = post_curation.get("source_taxon_db_xref")
    if source_taxon_db_xref and source_taxon_db_xref["source_taxon_db_xref_changed"]:
        before = source_taxon_db_xref["source_taxon_db_xrefs_before"]
        after = source_taxon_db_xref["source_taxon_db_xrefs_after"]
        before_text = ", ".join(before) if before else "none"
        after_text = ", ".join(after) if after else "none"
        lines.append(
            f"- source taxon db_xref: changed from {before_text} to {after_text}"
        )

    organism_taxonomy = post_curation.get("organism_taxonomy")
    if organism_taxonomy:
        before = organism_taxonomy["organism_taxonomy_before"]
        after = organism_taxonomy["organism_taxonomy_after"]
        taxid = organism_taxonomy["organism_taxonomy_taxid"]
        source = organism_taxonomy["organism_taxonomy_source"]
        source_text = source or "configured taxonomy lineage"
        taxid_text = f" for taxid {taxid}" if taxid else ""
        if organism_taxonomy["organism_taxonomy_changed"]:
            before_text = f'"{before}"' if before else "none"
            lines.append(
                "- ORGANISM taxonomy lineage: changed from "
                f"{before_text} to \"{after}\" ({source_text}{taxid_text})"
            )
        elif after and not organism_taxonomy["organism_taxonomy_skipped_reason"]:
            lines.append(
                "- ORGANISM taxonomy lineage: already matched "
                f"{source_text}{taxid_text}"
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

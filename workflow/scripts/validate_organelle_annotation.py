#!/usr/bin/env python3
"""Read-only FASTA/GenBank and CDS-translation consistency validation."""

import argparse
import csv
import hashlib
import json
import re
import warnings
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path

try:
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
    from Bio.SeqFeature import ExactPosition, Location
except ImportError:  # Checked by require_biopython at runtime.
    CodonTable = Seq = ExactPosition = Location = None

from reference_cds_qc import (
    SeqIO, feature_location_string, normalized_genbank_text_for_biopython,
    record_locus_id, require_biopython,
)


ISSUE_FIELDS = ["severity", "record_id", "feature_index", "feature_type", "gene",
                "location", "code", "message"]
AMINO_ACIDS = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))
AMINO_ACIDS.update(SEC="U", PYL="O", TERM="*", OTHER="X", XAA="X", ASX="B", GLX="Z")


def first(qualifiers, key, default=""):
    return (qualifiers.get(key) or [default])[0]


def location_bases(location, record_id):
    """Coordinates in biological extraction order, retaining remote record IDs."""
    return [(part.ref or record_id, position, part.strand)
            for part in location.parts for position in part]


def canonical_location(text, length):
    location = Location.fromstring(text, length=length)
    return feature_location_string(type("Feature", (), {"location": location})(), length)


def editing_index(path, records, emit):
    if path is None:
        return {}, {}
    try:
        document = json.loads(path.read_text())
        calls = document["feature_decisions"]
        if not isinstance(calls, list):
            raise ValueError("feature_decisions must be a list")
    except (OSError, ValueError, KeyError, TypeError) as error:
        emit("error", "editing_decisions_unreadable", str(error))
        return {}, {}
    indexed = defaultdict(list)
    for call in calls:
        try:
            record_id = call["record"]
            if record_id not in records:
                continue  # Removed features/records need not survive manual curation.
            key = (record_id, call.get("gene", "").casefold(),
                   canonical_location(call["location"], len(records[record_id])))
            indexed[key].append(call)
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            emit("error", "invalid_editing_decision", str(error))
    return indexed, document


def replay_editing(sequence, positions, calls, records, emit):
    applied = [c.get("accepted_sites", []) for c in calls if c.get("translation_added")]
    if not applied:
        return sequence, 0, False
    if any(not isinstance(sites, list) for sites in applied):
        emit("error", "invalid_editing_decision", "accepted_sites must be a list")
        return sequence, 0, True
    if any(sites != applied[0] for sites in applied[1:]):
        emit("error", "ambiguous_editing_decisions", "Conflicting decisions for the same CDS location")
        return sequence, 0, True
    edited = list(sequence)
    seen = set()
    for site in applied[0]:
        try:
            index = int(site["cds_pos"]) - 1
            if not 0 <= index < len(positions) or index in seen:
                raise ValueError("CDS position is out of range or repeated")
            seen.add(index)
            record_id, position, strand = positions[index]
            original = str(records[record_id].seq[position]).upper()
            expected = "G" if strand == -1 else "C"
            replacement = "A" if strand == -1 else "T"
            if (position + 1 != int(site["genomic_pos"])
                    or original != site["genomic_base"]
                    or original != expected or site["edited_base"] != replacement
                    or sequence[index] != "C"):
                raise ValueError("C-to-U call does not match the CDS/genomic coordinate, strand, or base")
            edited[index] = "T"
        except (KeyError, TypeError, ValueError, IndexError) as error:
            emit("error", "invalid_editing_site", str(error))
    return "".join(edited), sum(a != b for a, b in zip(sequence, edited)), True


def translation_exceptions(values, positions, coding, table, record_id, record_length, emit):
    replacements = {}
    partial_stop = False
    supported = True
    coordinates = [(r, p) for r, p, _strand in positions]
    for value in values:
        match = re.fullmatch(r"\(\s*pos:(.+),\s*aa:([A-Za-z*]+)\s*\)", value)
        if not match:
            emit("error", "invalid_transl_except", f"Invalid transl_except: {value}")
            supported = False
            continue
        location_text, amino_acid = match.groups()
        amino_acid = AMINO_ACIDS.get(amino_acid.upper())
        if amino_acid is None:
            emit("warning", "unsupported_transl_except", f"Unknown amino acid in {value}")
            supported = False
            continue
        try:
            location = Location.fromstring(location_text, length=record_length)
            sites = [(r, p) for r, p, _strand in location_bases(location, record_id)]
            wanted = set(sites)
            candidates = [i for i in range(0, len(coordinates), 3)
                          if len(coordinates[i:i+3]) == len(sites)
                          and set(coordinates[i:i+3]) == wanted]
            if len(candidates) != 1 or not 1 <= len(sites) <= 3:
                raise ValueError("Exception must identify one codon in the CDS reading frame")
            offset = candidates[0]
            if len(sites) < 3:
                if (amino_acid != "*" or offset + len(sites) != len(coding)
                        or coding[offset:] + "A" * (3-len(sites)) not in table.stop_codons):
                    raise ValueError("Partial TERM must be a terminal stop completed by polyadenylation")
                partial_stop = True
            codon_index = offset // 3
            if codon_index in replacements and replacements[codon_index] != amino_acid:
                raise ValueError("Conflicting exceptions for the same codon")
            replacements[codon_index] = amino_acid
        except (ValueError, TypeError) as error:
            emit("error", "invalid_transl_except", f"{value}: {error}")
            supported = False
    return replacements, partial_stop, supported


def validate_cds(feature, record, records, calls, default_table, emit):
    q = feature.qualifiers
    result = {"translation_status": "not_checked", "editing_sites_replayed": 0}
    if "pseudo" in q or "pseudogene" in q:
        result["translation_status"] = "skipped_pseudogene"
        return result
    try:
        offset = int(first(q, "codon_start", "1")) - 1
        table_id = int(first(q, "transl_table", default_table))
        if offset not in (0, 1, 2):
            raise ValueError("codon_start must be 1, 2, or 3")
        table = CodonTable.unambiguous_dna_by_id[table_id]
    except (ValueError, TypeError, KeyError) as error:
        emit("error", "invalid_translation_parameters", str(error))
        return result
    if getattr(feature.location, "operator", "join") != "join":
        emit("warning", "unsupported_location_operator", "Only join CDS locations can be reconstructed")
        return result
    record_id = record_locus_id(record)
    try:
        sequence = str(feature.extract(record.seq, references={k: r.seq for k, r in records.items()})).upper()
        positions = location_bases(feature.location, record_id)
    except (ValueError, TypeError, KeyError) as error:
        emit("error", "CDS_extraction_failed", str(error))
        return result
    sequence, replayed, has_decision = replay_editing(sequence, positions, calls, records, emit)
    result.update(editing_sites_replayed=replayed, transl_table=table_id, codon_start=offset+1,
                  extracted_length=len(sequence))
    coding = sequence[offset:]
    first_part, last_part = feature.location.parts[0], feature.location.parts[-1]
    partial5 = not isinstance(first_part.end if first_part.strand == -1 else first_part.start, ExactPosition)
    partial3 = not isinstance(last_part.start if last_part.strand == -1 else last_part.end, ExactPosition)
    result.update(partial5=partial5, partial3=partial3)
    exceptions = {str(v).casefold() for v in q.get("exception", [])}
    unmodeled = exceptions - {"rna editing", "trans-splicing", "trans_splicing"}
    if "rna editing" in exceptions and not has_decision:
        unmodeled.add("RNA editing without matching applied-site decisions")
    if "ribosomal_slippage" in q and len(feature.location.parts) == 1:
        unmodeled.add("ribosomal slippage without a reconstructable join")
    for exception in sorted(unmodeled):
        emit("warning", "unverified_translation_exception", exception)
    replacements, partial_stop, supported = translation_exceptions(
        q.get("transl_except", []), positions[offset:], coding, table,
        record_id, len(record), emit)
    stored = "".join(str(first(q, "translation")).split())
    if stored.endswith("*"):
        stored = stored[:-1]
    if not stored:
        emit("warning", "translation_missing", "No stored translation to compare")
    if not supported or unmodeled:
        result["translation_status"] = "review_exception"
        return result
    ambiguous = bool(re.search(r"[^ACGT]", coding))
    if ambiguous:
        emit("warning", "ambiguous_CDS_bases", "CDS contains ambiguous/non-ACGT bases")
    problems = []
    if len(coding) < 3:
        problems.append(("CDS_too_short", "No complete codon after codon_start"))
    if len(coding) % 3 and not partial3 and not partial_stop:
        problems.append(("CDS_length_not_in_frame", "Complete 3-prime CDS length is not divisible by three"))
    try:
        amino_acids = list(str(Seq(coding[:len(coding)//3*3]).translate(table=table_id)))
    except (ValueError, CodonTable.TranslationError) as error:
        emit("error", "CDS_translation_failed", str(error))
        return result
    if not partial5 and coding[:3] in table.start_codons and amino_acids:
        amino_acids[0] = "M"
    elif not partial5 and 0 not in replacements and not ambiguous:
        problems.append(("invalid_start_codon", "Complete 5-prime CDS lacks a valid initiation codon"))
    # Some genetic codes have context-dependent stop codons; CDS termination takes priority.
    if not partial3 and not len(coding) % 3 and coding[-3:] in table.stop_codons and amino_acids:
        amino_acids[-1] = "*"
    for index, amino_acid in replacements.items():
        if index == len(amino_acids) and partial_stop:
            amino_acids.append(amino_acid)
        elif index < len(amino_acids):
            amino_acids[index] = amino_acid
    if not partial3 and (not amino_acids or amino_acids[-1] != "*") and not ambiguous:
        problems.append(("terminal_stop_missing", "Complete 3-prime CDS lacks a terminal stop"))
    if amino_acids and amino_acids[-1] == "*":
        amino_acids.pop()
    expected = "".join(amino_acids)
    if "*" in expected:
        problems.append(("internal_stop", "Reconstructed protein contains an internal stop"))
    if stored and stored != expected:
        problems.append(("translation_mismatch", "Stored translation differs from reconstructed CDS"))
    for code, message in problems:
        emit("error", code, message)
    result.update(translation_status="invalid" if problems else "matched" if stored else "missing",
                  reconstructed_protein_length=len(expected), stored_protein_length=len(stored))
    return result


def validate_annotation(*, annotation, fasta, decisions=None, default_transl_table="1"):
    require_biopython()
    issues, cds_rows = [], []

    def emit(severity, code, message, **context):
        issues.append({**dict.fromkeys(ISSUE_FIELDS, ""), **context,
                       "severity": severity, "code": code, "message": message})

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            records = list(SeqIO.parse(StringIO(normalized_genbank_text_for_biopython(annotation)), "genbank"))
        fastas = list(SeqIO.parse(fasta, "fasta"))
    except (OSError, ValueError, TypeError) as error:
        emit("error", "input_parse_error", str(error))
        records, fastas = [], []
    if not records or not fastas:
        emit("error", "empty_input", "Both FASTA and GenBank must contain records")
    record_ids = [record_locus_id(r) for r in records]
    fasta_ids = [r.id for r in fastas]
    for label, ids in [("GenBank", record_ids), ("FASTA", fasta_ids)]:
        for identifier, count in Counter(ids).items():
            if count > 1:
                emit("error", "duplicate_record_id", f"{label}: {identifier} occurs {count} times")
    if set(record_ids) != set(fasta_ids):
        emit("error", "record_id_mismatch", f"GenBank-only: {sorted(set(record_ids)-set(fasta_ids))}; FASTA-only: {sorted(set(fasta_ids)-set(record_ids))}")
    by_id = dict(zip(record_ids, records))
    fasta_by_id = {r.id: r for r in fastas}
    indexed, document = editing_index(decisions, by_id, emit)
    table_default = document.get("transl_table_default", default_transl_table)
    for record in records:
        record_id = record_locus_id(record)
        if record_id in fasta_by_id and str(record.seq).upper() != str(fasta_by_id[record_id].seq).upper():
            emit("error", "sequence_mismatch", "GenBank bases differ from the annotation FASTA", record_id=record_id)
        for index, feature in enumerate(record.features):
            context = dict(record_id=record_id, feature_index=index, feature_type=feature.type,
                           gene=str(first(feature.qualifiers, "gene")))

            def feature_issue(severity, code, message):
                emit(severity, code, message, **context)

            location = feature.location
            if location is None:
                feature_issue("error", "unparsed_location", "Feature location could not be parsed")
                if feature.type == "CDS":
                    cds_rows.append({**context, "translation_status": "not_checked"})
                continue
            try:
                context["location"] = feature_location_string(feature, len(record))
                valid = True
                for part in location.parts:
                    source = by_id.get(part.ref or record_id)
                    if source is None:
                        feature_issue("warning", "unresolved_remote_location", f"Record not provided: {part.ref}")
                        valid = False
                    elif not 0 <= int(part.start) <= int(part.end) <= len(source):
                        feature_issue("error", "location_out_of_bounds", "Feature part is outside its sequence")
                        valid = False
                    elif feature.type in {"CDS", "tRNA", "rRNA"} and int(part.start) == int(part.end):
                        feature_issue("error", "empty_product_location", "Product feature has an empty location")
                        valid = False
            except (ValueError, TypeError) as error:
                feature_issue("error", "invalid_location", str(error))
                valid = False
            if feature.type == "CDS":
                key = (record_id, context["gene"].casefold(), context.get("location", ""))
                details = validate_cds(feature, record, by_id, indexed.get(key, []), table_default, feature_issue) if valid else {"translation_status": "not_checked"}
                cds_rows.append({**context, **details})
    counts = Counter(i["severity"] for i in issues)
    return {"schema_version": 1, "annotation": str(annotation), "fasta": str(fasta),
            "rna_editing_decisions": str(decisions) if decisions else None,
            "input_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in (annotation, fasta, decisions) if p is not None and p.is_file()},
            "record_count": len(records), "CDS_count": len(cds_rows),
            "error_count": counts["error"], "warning_count": counts["warning"],
            "status": "error" if counts["error"] else "review" if counts["warning"] else "pass",
            "issues": issues, "CDS": cds_rows}


def run_validation(*, annotation, fasta, report_json, report_tsv, decisions=None, default_transl_table="1"):
    report = validate_annotation(annotation=annotation, fasta=fasta, decisions=decisions,
                                 default_transl_table=default_transl_table)
    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report_tsv.parent.mkdir(parents=True, exist_ok=True)
    with report_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ISSUE_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(report["issues"])
    return report


def check_report(path):
    report = json.loads(path.read_text())
    if report["error_count"]:
        raise ValueError(f"Organelle annotation has {report['error_count']} validation error(s); see {path}")
    for name, expected in report["input_sha256"].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Validation report is stale; input changed: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-report", type=Path, help="Reject error or stale reports before release")
    parser.add_argument("--annotation", type=Path)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--decisions", type=Path)
    parser.add_argument("--default-transl-table", default="1")
    parser.add_argument("--report-json", type=Path)
    parser.add_argument("--report-tsv", type=Path)
    args = parser.parse_args()
    if args.check_report:
        try:
            check_report(args.check_report)
        except (OSError, ValueError, KeyError) as error:
            parser.exit(1, f"{error}\n")
        return
    if any(getattr(args, k) is None for k in ("annotation", "fasta", "report_json", "report_tsv")):
        parser.error("--annotation, --fasta, --report-json and --report-tsv are required")
    values = vars(args)
    values.pop("check_report")
    report = run_validation(**values)
    print(f"Validation: {report['status']}; {report['error_count']} error(s), {report['warning_count']} review warning(s)")


if __name__ == "__main__":
    main()

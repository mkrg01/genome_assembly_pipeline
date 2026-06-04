import argparse
import json
from pathlib import Path

from organelle_annotation_utils import (
    normalize_taxid,
    taxonomy_record_from_ncbi_taxa,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare a GenBank-ready NCBI taxonomy lineage record."
    )
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--lineage", type=Path, required=True)
    return parser.parse_args()


def ensure_ete_cache_dir():
    cache_dir = Path.home() / ".local" / "share" / "ete"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def main():
    args = parse_args()
    taxid = normalize_taxid(args.taxid)

    args.db.parent.mkdir(parents=True, exist_ok=True)
    args.lineage.parent.mkdir(parents=True, exist_ok=True)
    ensure_ete_cache_dir()

    from ete4 import NCBITaxa

    ncbi = NCBITaxa(
        dbfile=str(args.db),
        update=True,
    )
    record = taxonomy_record_from_ncbi_taxa(ncbi, taxid)
    record.update(
        {
            "taxonomy_database": str(args.db),
        }
    )
    args.lineage.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

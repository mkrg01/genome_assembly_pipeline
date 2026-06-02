import json
import shutil
from pathlib import Path


GENBANK_SUFFIXES = (".gb", ".gbk", ".gbf")


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


def write_run_manifest(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")

"""Download and verify an immutable BRAKER4 source archive."""

import argparse
import hashlib
import json
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def install_archive(archive, source, expected_sha256):
    actual = sha256_file(archive)
    if actual != expected_sha256:
        raise ValueError(f"BRAKER4 SHA256 mismatch: expected {expected_sha256}, got {actual}")
    source = Path(source)
    source.parent.mkdir(parents=True, exist_ok=True)
    # Source directories are version-specific. Do not silently replace a local edit.
    if source.exists():
        raise FileExistsError(f"BRAKER4 source destination already exists: {source}")
    with tempfile.TemporaryDirectory(dir=source.parent) as tmp:
        unpacked = Path(tmp)
        with tarfile.open(archive) as bundle:
            bundle.extractall(unpacked, filter="data")
        roots = list(unpacked.iterdir())
        if len(roots) != 1 or not (roots[0] / "Snakefile").is_file():
            raise ValueError("BRAKER4 archive must contain one source directory with a Snakefile")
        roots[0].rename(source)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    args.archive.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.archive.parent) as tmp:
        downloaded = Path(tmp) / "source.tar.gz"
        with urllib.request.urlopen(args.url, timeout=120) as response, downloaded.open("wb") as out:
            shutil.copyfileobj(response, out)
        install_archive(downloaded, args.source, args.sha256)
        downloaded.replace(args.archive)
    args.manifest.write_text(json.dumps({
        "tool": "BRAKER4", "version": args.version, "commit": args.commit,
        "url": args.url, "sha256": args.sha256,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()

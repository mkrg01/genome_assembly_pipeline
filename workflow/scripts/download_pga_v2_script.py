import argparse
import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Download a pinned PGA v2.0 script.")
    parser.add_argument("--url", required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    return parser.parse_args()


def sha256sum(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=output_path.name + ".", suffix=".tmp", dir=output_path.parent, delete=False
    ) as tmp_handle:
        tmp_path = Path(tmp_handle.name)
        with urllib.request.urlopen(url) as response:
            shutil.copyfileobj(response, tmp_handle)
    tmp_path.replace(output_path)


def write_manifest(args, actual_sha256):
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(
            {
                "tool": "pga_v2",
                "source": "github",
                "commit": args.commit,
                "url": args.url,
                "expected_sha256": args.sha256,
                "actual_sha256": actual_sha256,
                "script": args.output.as_posix(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def main():
    args = parse_args()
    if not args.output.exists() or sha256sum(args.output) != args.sha256:
        download_file(args.url, args.output)

    actual_sha256 = sha256sum(args.output)
    if actual_sha256 != args.sha256:
        raise RuntimeError(
            f"SHA256 mismatch for {args.output}. Expected {args.sha256}, got {actual_sha256}."
        )

    current_mode = args.output.stat().st_mode
    os.chmod(args.output, current_mode | 0o755)
    write_manifest(args, actual_sha256)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
import sys
import hashlib
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Validate MD5 checksum of a file against a provided MD5 file.")
    parser.add_argument("file", type=Path, required=True, help="Path to the file to validate.")
    parser.add_argument("md5_file", type=Path, required=True, help="Path to the MD5 file containing expected checksum.")
    return parser.parse_args()

def check_md5(file_path, md5_path):
    with open(md5_path) as f:
        expected_md5, expected_name = f.readline().strip().split()
        expected_md5 = expected_md5.lower()

    with open(file_path, "rb") as f:
        actual_md5 = hashlib.md5(f.read()).hexdigest().lower()

    if actual_md5 != expected_md5:
        raise ValueError(f"MD5 mismatch!\nExpected: {expected_md5}\nActual: {actual_md5}")
    else:
        print("MD5 check passed.")

if __name__ == "__main__":
    args = parse_args()

    file_path = args.file
    md5_path = args.md5_file

    if not file_path.exists() or not md5_path.exists():
        print("Error: One or both files do not exist.", file=sys.stderr)
        sys.exit(1)

    try:
        check_md5(file_path, md5_path)
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)

import argparse
import gzip
import shutil
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Copy a plain-text or gzipped submission file into a normalized gzip-compressed output."
        )
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def open_input(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def copy_to_gzip(input_path: Path, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_input(input_path) as in_handle, gzip.open(output_path, "wb") as out_handle:
        shutil.copyfileobj(in_handle, out_handle)


def main():
    args = parse_args()
    copy_to_gzip(args.input, args.output)


if __name__ == "__main__":
    main()

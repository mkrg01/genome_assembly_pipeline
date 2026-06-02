import argparse
import hashlib
import json
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path


FIGSHARE_FILES_API = "https://api.figshare.com/v2/articles/{article_id}/versions/{version}/files"


def parse_args():
    parser = argparse.ArgumentParser(description="Download and unpack the pinned PMGA bundle.")
    parser.add_argument("--article-id", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--file-name", required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def fetch_json(url):
    with urllib.request.urlopen(url) as response:
        return json.loads(response.read().decode("utf-8"))


def select_file(files, file_name):
    for item in files:
        if item.get("name") == file_name:
            return item
    names = ", ".join(sorted(item.get("name", "<unnamed>") for item in files))
    raise RuntimeError(f"Could not find '{file_name}' in Figshare files: {names}")


def md5sum(path):
    digest = hashlib.md5()
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


def top_level_members(tar):
    names = set()
    for member in tar.getmembers():
        first = Path(member.name).parts[0] if Path(member.name).parts else None
        if first and first not in (".", ".."):
            names.add(first)
    return sorted(names)


def checked_members(tar, destination):
    destination = destination.resolve()
    for member in tar.getmembers():
        if member.issym() or member.islnk():
            raise RuntimeError(f"Unsafe link in PMGA archive: {member.name}")
        member_path = (destination / member.name).resolve()
        if destination not in [member_path, *member_path.parents]:
            raise RuntimeError(f"Unsafe path in PMGA archive: {member.name}")
        yield member


def unpack_archive(archive, bundle):
    bundle_parent = bundle.parent
    bundle_parent.mkdir(parents=True, exist_ok=True)
    if bundle.exists():
        shutil.rmtree(bundle)

    with tempfile.TemporaryDirectory(dir=bundle_parent) as tmpdir:
        tmpdir = Path(tmpdir)
        with tarfile.open(archive) as tar:
            members = top_level_members(tar)
            tar.extractall(tmpdir, members=checked_members(tar, tmpdir))

        preferred = tmpdir / bundle.name
        if preferred.exists():
            shutil.move(str(preferred), bundle)
        elif len(members) == 1 and (tmpdir / members[0]).exists():
            shutil.move(str(tmpdir / members[0]), bundle)
        else:
            bundle.mkdir(parents=True, exist_ok=True)
            for child in tmpdir.iterdir():
                shutil.move(str(child), bundle / child.name)


def write_manifest(args, file_info, actual_md5):
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(
            {
                "tool": "pmga",
                "source": "figshare",
                "article_id": args.article_id,
                "version": args.version,
                "file_name": args.file_name,
                "file_id": file_info.get("id"),
                "download_url": file_info.get("download_url"),
                "size": file_info.get("size"),
                "expected_md5": file_info.get("computed_md5") or file_info.get("supplied_md5"),
                "actual_md5": actual_md5,
                "archive": args.archive.as_posix(),
                "bundle": args.bundle.as_posix(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def main():
    args = parse_args()
    url = FIGSHARE_FILES_API.format(article_id=args.article_id, version=args.version)
    file_info = select_file(fetch_json(url), args.file_name)
    expected_md5 = file_info.get("computed_md5") or file_info.get("supplied_md5")
    if not expected_md5:
        raise RuntimeError(f"Figshare metadata for {args.file_name} does not include an MD5.")

    archive_was_downloaded = False
    if not args.archive.exists() or md5sum(args.archive) != expected_md5:
        download_file(file_info["download_url"], args.archive)
        archive_was_downloaded = True

    actual_md5 = md5sum(args.archive)
    if actual_md5 != expected_md5:
        raise RuntimeError(
            f"MD5 mismatch for {args.archive}. Expected {expected_md5}, got {actual_md5}."
        )

    if archive_was_downloaded or not args.bundle.exists():
        unpack_archive(args.archive, args.bundle)
    write_manifest(args, file_info, actual_md5)


if __name__ == "__main__":
    main()

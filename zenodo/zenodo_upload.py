#!/usr/bin/env python3
"""
zenodo_upload.py — upload the tarballs in zenodo/archives/ to a new Zenodo deposition.

Prerequisites:
    1. Run `bash zenodo/make_archives.sh` first (creates zenodo/archives/*.tar.gz).
    2. Create a Zenodo personal access token (Zenodo > Applications > Personal access
       tokens) with the `deposit:write` and `deposit:actions` scopes.
    3. export ZENODO_TOKEN=xxxxxxxx

Usage:
    python zenodo/zenodo_upload.py                 # upload to the real Zenodo
    python zenodo/zenodo_upload.py --sandbox       # upload to sandbox.zenodo.org (test)
    python zenodo/zenodo_upload.py --publish       # also publish (mint the DOI) at the end

This creates a draft deposition and uploads every archive. It prints the reserved DOI.
Review the draft on the Zenodo website, add metadata (authors, title, description),
then publish there — or pass --publish to publish immediately.
"""
import argparse
import os
import sys
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("Please `pip install requests` first.")

ARCHIVE_DIR = Path(__file__).resolve().parent / "archives"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org")
    ap.add_argument("--publish", action="store_true", help="publish (mint DOI) after upload")
    args = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("Set ZENODO_TOKEN environment variable first.")

    base = "https://sandbox.zenodo.org" if args.sandbox else "https://zenodo.org"
    params = {"access_token": token}

    tarballs = sorted(ARCHIVE_DIR.glob("*.tar.gz"))
    if not tarballs:
        sys.exit(f"No archives in {ARCHIVE_DIR}. Run make_archives.sh first.")

    # 1. create an empty deposition
    r = requests.post(f"{base}/api/deposit/depositions", params=params, json={})
    r.raise_for_status()
    dep = r.json()
    dep_id = dep["id"]
    bucket = dep["links"]["bucket"]
    doi = dep.get("metadata", {}).get("prereserve_doi", {}).get("doi", "(reserved on save)")
    print(f"Created deposition {dep_id}  reserved DOI: {doi}")
    print(f"Draft URL: {dep['links'].get('html')}")

    # 2. stream each tarball into the bucket
    for tb in tarballs:
        print(f"  uploading {tb.name} ({tb.stat().st_size/1e9:.2f} GB) ...", flush=True)
        with open(tb, "rb") as fh:
            up = requests.put(f"{bucket}/{tb.name}", data=fh, params=params)
        up.raise_for_status()
    # upload the checksum file too
    sums = ARCHIVE_DIR / "SHA256SUMS"
    if sums.exists():
        with open(sums, "rb") as fh:
            requests.put(f"{bucket}/SHA256SUMS", data=fh, params=params).raise_for_status()

    print(f"\nAll files uploaded to deposition {dep_id}.")
    if args.publish:
        pub = requests.post(f"{base}/api/deposit/depositions/{dep_id}/actions/publish",
                            params=params)
        pub.raise_for_status()
        print("Published. Final DOI:", pub.json()["doi"])
    else:
        print("Draft NOT published. Add metadata and publish at:", dep["links"].get("html"))
    print("\nAfter publishing, paste the DOI into README.md and DATA.md.")


if __name__ == "__main__":
    main()

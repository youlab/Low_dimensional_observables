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
    python zenodo/zenodo_upload.py --deposition ID # resume: upload into an existing draft,
                                                   #   skipping files already present

This creates a draft deposition and uploads every archive. It prints the reserved DOI.
Review the draft on the Zenodo website, add metadata (authors, title, description),
then publish there — or pass --publish to publish immediately.

The upload is resumable and retries transient server errors (500/502/503/504) and
dropped connections with exponential backoff, so a large multi-GB transfer survives the
occasional Zenodo hiccup. If a run dies anyway, re-run with --deposition <ID> (printed at
the top of every run) to continue into the same draft without re-sending finished files.
"""
import argparse
import os
import sys
import time
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("Please `pip install requests` first.")

ARCHIVE_DIR = Path(__file__).resolve().parent / "archives"

RETRY_STATUS = {500, 502, 503, 504}
MAX_ATTEMPTS = 10          # ride out sustained flaky-network windows
MAX_BACKOFF = 120          # seconds


def put_with_retry(url, path, params):
    """PUT a file, retrying transient server/connection errors with backoff."""
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with open(path, "rb") as fh:
                resp = requests.put(url, data=fh, params=params)
            if resp.status_code in RETRY_STATUS:
                raise requests.exceptions.HTTPError(
                    f"{resp.status_code} {resp.reason}", response=resp)
            resp.raise_for_status()
            return resp
        except (requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.HTTPError) as e:
            status = getattr(getattr(e, "response", None), "status_code", None)
            # non-retryable HTTP errors (e.g. 400/401/403) should fail immediately
            if status is not None and status not in RETRY_STATUS:
                raise
            if attempt == MAX_ATTEMPTS:
                raise
            wait = min(MAX_BACKOFF, 2 ** attempt)
            print(f"    transient error ({e}); retry {attempt}/{MAX_ATTEMPTS - 1} "
                  f"in {wait}s ...", flush=True)
            time.sleep(wait)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org")
    ap.add_argument("--publish", action="store_true", help="publish (mint DOI) after upload")
    ap.add_argument("--deposition", type=int, default=None,
                    help="resume into an existing draft deposition ID (skip files present)")
    args = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("Set ZENODO_TOKEN environment variable first.")

    base = "https://sandbox.zenodo.org" if args.sandbox else "https://zenodo.org"
    params = {"access_token": token}

    tarballs = sorted(ARCHIVE_DIR.glob("*.tar.gz"))
    if not tarballs:
        sys.exit(f"No archives in {ARCHIVE_DIR}. Run make_archives.sh first.")

    # 1. create a new deposition, or resume into an existing draft
    if args.deposition is not None:
        r = requests.get(f"{base}/api/deposit/depositions/{args.deposition}", params=params)
        r.raise_for_status()
        dep = r.json()
        if dep.get("submitted"):
            sys.exit(f"Deposition {args.deposition} is already published; cannot add files.")
        print(f"Resuming into existing deposition {dep['id']}")
    else:
        r = requests.post(f"{base}/api/deposit/depositions", params=params, json={})
        r.raise_for_status()
        dep = r.json()
    dep_id = dep["id"]
    bucket = dep["links"]["bucket"]
    doi = dep.get("metadata", {}).get("prereserve_doi", {}).get("doi", "(reserved on save)")
    already = {f["filename"] for f in dep.get("files", [])}
    print(f"Deposition {dep_id}  reserved DOI: {doi}")
    print(f"Draft URL: {dep['links'].get('html')}")
    if already:
        print(f"Already present ({len(already)}), will skip: {', '.join(sorted(already))}")

    # 2. stream each file into the bucket, skipping any already uploaded.
    #    A file that fails all retries is recorded but does NOT abort the run, so one
    #    stubborn upload can't block the others queued behind it. Re-run with
    #    --deposition <id> to retry only what is still missing.
    uploads = list(tarballs) + [ARCHIVE_DIR / n for n in ("SHA256SUMS", "README.md")]
    failed = []
    for item in uploads:
        if not item.exists():
            continue
        if item.name in already:
            print(f"  skipping {item.name} (already uploaded)", flush=True)
            continue
        print(f"  uploading {item.name} ({item.stat().st_size/1e9:.2f} GB) ...", flush=True)
        try:
            put_with_retry(f"{bucket}/{item.name}", item, params)
        except Exception as e:  # noqa: BLE001 — keep going, report at the end
            print(f"  FAILED {item.name}: {e}", flush=True)
            failed.append(item.name)

    if failed:
        print(f"\n{len(failed)} file(s) FAILED after retries: {', '.join(failed)}")
        print(f"Re-run to retry only these:  "
              f"sbatch zenodo/up_zenodo.sh --deposition {dep_id}")
        sys.exit(1)

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

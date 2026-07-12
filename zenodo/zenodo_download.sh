#!/usr/bin/env bash
# zenodo_download.sh — fetch the archived artifacts from Zenodo, verify checksums, and
# extract them back into their original folders so the notebooks can recompute results.
#
# Usage (from the repo root):
#     bash zenodo/zenodo_download.sh 10.5281/zenodo.XXXXXXX          # all tarballs
#     bash zenodo/zenodo_download.sh 10.5281/zenodo.XXXXXXX glv_models glv_saved_sims
#
# Pass a record DOI (or a bare record id). With no extra args it downloads every file
# in the record; otherwise only the named tarballs (without the .tar.gz suffix).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

DOI="${1:-}"
[ -z "$DOI" ] && { echo "Usage: $0 <zenodo-doi-or-record-id> [tarball names...]"; exit 1; }
shift || true
RECORD_ID="${DOI##*.}"   # 10.5281/zenodo.123 -> 123 ; also works for a bare id
DEST="zenodo/archives"
mkdir -p "$DEST"

echo "Querying Zenodo record $RECORD_ID ..."
API="https://zenodo.org/api/records/$RECORD_ID"
# list files: name + download url
mapfile -t FILES < <(python3 - "$API" "$@" <<'PY'
import json, sys, urllib.request
api, *want = sys.argv[1], sys.argv[2:]
data = json.load(urllib.request.urlopen(sys.argv[1]))
for f in data.get("files", []):
    name = f.get("key") or f.get("filename")
    url  = f["links"].get("self") or f["links"].get("download")
    if want and name.removesuffix(".tar.gz") not in want and name != "SHA256SUMS":
        continue
    print(f"{name}\t{url}")
PY
)

for line in "${FILES[@]}"; do
  name="${line%%$'\t'*}"; url="${line##*$'\t'}"
  echo "  downloading $name ..."
  curl -fSL "$url" -o "$DEST/$name"
done

if [ -f "$DEST/SHA256SUMS" ]; then
  echo "Verifying checksums ..."
  ( cd "$DEST" && sha256sum -c --ignore-missing SHA256SUMS )
fi

echo "Extracting into repo (paths are repo-relative) ..."
for tb in "$DEST"/*.tar.gz; do
  [ -e "$tb" ] || continue
  echo "  tar xzf $tb"
  tar xzf "$tb" -C "$REPO_ROOT"
done
echo "Done. Artifacts restored to their original folders."

#!/usr/bin/env bash
# zenodo_download.sh — fetch the archived artifacts from Zenodo, verify checksums, and
# extract them back into their original folders so the notebooks can recompute results.
#
# Usage (from the repo root):
#     bash zenodo/zenodo_download.sh 10.5281/zenodo.21368291          # all tarballs
#     bash zenodo/zenodo_download.sh 10.5281/zenodo.21368291 glv_models glv_saved_sims
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
# Read into an array with a while-loop (portable; macOS bash 3.2 has no `mapfile`).
FILES=()
while IFS= read -r line; do
  [ -n "$line" ] && FILES+=("$line")
done < <(python3 - "$API" "$@" <<'PY'
import json, sys, urllib.request
api, want = sys.argv[1], sys.argv[2:]
data = json.load(urllib.request.urlopen(sys.argv[1]))
for f in data.get("files", []):
    name = f.get("key") or f.get("filename")
    url  = f["links"].get("self") or f["links"].get("download")
    if want and name.removesuffix(".tar.gz") not in want and name != "SHA256SUMS":
        continue
    print(f"{name}\t{url}")
PY
)

[ "${#FILES[@]}" -eq 0 ] && { echo "No matching files found in record $RECORD_ID."; exit 1; }
for line in "${FILES[@]}"; do
  name="${line%%$'\t'*}"; url="${line##*$'\t'}"
  echo "  downloading $name ..."
  curl -fSL "$url" -o "$DEST/$name"
done

if [ -f "$DEST/SHA256SUMS" ]; then
  echo "Verifying checksums ..."
  if command -v sha256sum >/dev/null 2>&1; then
    ( cd "$DEST" && sha256sum -c --ignore-missing SHA256SUMS )
  else
    # macOS has no sha256sum; shasum -a 256 lacks --ignore-missing, so feed it only present files.
    ( cd "$DEST" \
      && while read -r sum name; do [ -f "$name" ] && printf '%s  %s\n' "$sum" "$name"; done < SHA256SUMS \
      | shasum -a 256 -c - )
  fi
fi

echo "Extracting into repo (paths are repo-relative) ..."
for tb in "$DEST"/*.tar.gz; do
  [ -e "$tb" ] || continue
  echo "  tar xzf $tb"
  tar xzf "$tb" -C "$REPO_ROOT"
done
echo "Done. Artifacts restored to their original folders."

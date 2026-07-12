# Zenodo archive helpers

Scripts to bundle, upload, and restore the heavy artifacts (trained models + simulation
arrays) that are kept out of git. See `../DATA.md` for the tarball → folder manifest.

## Upload (maintainer, once)
```bash
bash zenodo/make_archives.sh                 # -> zenodo/archives/*.tar.gz + SHA256SUMS
pip install requests                         # if needed
export ZENODO_TOKEN=xxxxxxxx                 # token with deposit:write, deposit:actions
python zenodo/zenodo_upload.py               # creates a draft; --sandbox to test first
# add authors/title/description on the Zenodo draft page, then Publish (or use --publish)
# finally: paste the DOI into README.md and DATA.md
```

## Download (anyone reproducing the work)
```bash
bash zenodo/zenodo_download.sh 10.5281/zenodo.XXXXXXX          # all artifacts
bash zenodo/zenodo_download.sh 10.5281/zenodo.XXXXXXX glv_models   # just one tarball
```

`zenodo/archives/` is git-ignored — it only holds the local tarballs during upload/download.

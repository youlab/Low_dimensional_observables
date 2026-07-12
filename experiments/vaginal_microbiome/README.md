# Vaginal microbiome (24 subjects)

Cross-validated VAE analysis of 24-subject vaginal microbiome time series; effective observable dimension Eᴄ per subject.

## Model
`VAE_for_vaginal_microbiome_CV.py`

## Run order
1. `VAE_for_vaginal_microbiome_CV.py` / `Train_vag_microbiome_CV.py` — train CV VAEs; `VAE_all_24subjects.py` / `Retrain_vag_microbiome_Ec.py` for full models (`retrain_vag_full.sh`).
2. `Vaginal_microbiome_Ec_CV.ipynb` — Eᴄ / FUV figures (reads `vaginal_embedding_FUV_CV/` caches).

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `vae_models_CV/`, `vaginal_dataset_CV/`
- **In this repo** (needed to reproduce the figures without a download): `vaginal_embedding_FUV_CV/`, `meta_data/`

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

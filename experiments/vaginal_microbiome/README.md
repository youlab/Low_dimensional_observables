# Vaginal microbiome (24 subjects)

Cross-validated VAE analysis of 24-subject vaginal microbiome time series; effective observable dimension Eᴄ.

## Model
`VAE_for_vaginal_microbiome_CV.py`

## Run order
1. `VAE_for_vaginal_microbiome_CV.py` / `Train_vag_microbiome_CV.py` — train CV VAEs (`sbatch train_vag_CV.sh <N_TARGET>` sweeps the 24 held-out subjects × 3 trials).
2. `Vaginal_microbiome_Ec_CV.ipynb` — Eᴄ / FUV figures (reads `vaginal_embedding_FUV_CV/` caches).

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `vae_models_CV/`, `vaginal_dataset_CV/`
- **In this repo** (needed to reproduce the figures without a download): `vaginal_embedding_FUV_CV/`
- **Provenance** (not read by any script): `meta_data/microbiome_timeseries_ranked.npz` — the per-subject relative-abundance series that `vaginal_dataset_CV/` was built from by sliding-window augmentation (14-day windows, step 1 day; taxa kept when their cumulative relative abundance over the window exceeds 0.1; up to 100 taxon combinations per window).

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

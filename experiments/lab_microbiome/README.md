# Lab microbiome (soil & water communities)

Cross-validated VAE analysis of measured soil and water communities: effective observable dimension Eᴄ, community-collapse prediction, community-shift analysis, and OLS forecasting.

## Model
`VAE_lab_microbiome_CV.py`

## Run order
1. `VAE_lab_microbiome_CV.py` / `Train_lab_microbiome_CV.py` — train cross-validated VAEs; `VAE_lab_microbiome_full.py` / `Retrain_lab_microbiome.py` for full-data models (`retrain_full_VAE.sh`).
2. `LabMicrobiome_Ec_CV.ipynb` — Eᴄ / FUV figures (reads `lab_microbiome_embedding_FUV_CV/` caches; recompute needs the CV models from Zenodo).
3. `MLP_collapse_prediction.py` (`train_mlp_collapse.sh`) + `MLP_commshift.ipynb` — collapse prediction & community shift.
4. `OLS_raw_forecast.ipynb` / `OLS_delta_forecast.ipynb` (+ `plot_*_allmedia.py`) — OLS forecasting; outputs in `target_forecast_OLS_regression_raw/`, `delta_forecast_OLS/`.

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `vae_models_CV/`, `mlp_models/` (collapse weights), `lab_microbiome_datasets_CV/`
- **In this repo** (needed to reproduce the figures without a download): `lab_microbiome_embedding_FUV_CV/`, `community_shift/`, `collapse_forecast_results/`, `delta_forecast_OLS/`, `target_forecast_OLS_regression_raw/`, collapse `*_traindoc.pth`

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

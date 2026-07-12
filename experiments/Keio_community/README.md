# Engineered target populations (fluorescent proteins)

VAE analysis of measured growth dynamics of engineered fluorescent-protein (FP) populations and community-assembly experiments.

## Model
`VAE_model.py`, `MLPVAE.py`

## Run order
1. `train_exp_target.py` / `train_growth_curve_inference.py` (+ `MLPVAE.py`) — train VAEs per FP target (`Train_EGFP_curve.sh`, `Train_mCherry_curve.sh`, `Train_mTagBFP2_curve.sh`, `Train_LSSmOrange_curve.sh`, `Train_comm_assembly.sh`).
2. `VAE_on_FP_populations.ipynb` / `VAE_on_single_target.ipynb` — Eᴄ / reconstruction figures. Use the in-repo `Nano_target_dynamics_{train,test}.npy` splits.

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `vae_models/`, `mlp_models/`
- **In this repo** (needed to reproduce the figures without a download): `Nano_target_dynamics*.npy` (measured train/test splits), `sequenced_data/`

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

# Stochastic logistic growth

Independent logistic growers with intrinsic noise. Probes Eᴄ when dynamics are stochastic rather than deterministic.

## Model
`VAE_model.py`

## Run order
1. `simulation.ipynb` — generate trajectories (params saved as `saved_sims/I*/params_I*.pkl`).
2. `train_slogistic.py` — train VAEs (`VAE_slogistic_I1.sh`, `VAE_slogistic_I2-5.sh`, `VAE_slogistic_I6-9.sh`).
3. `Summary_stochastic_logistic.ipynb` — FUV / Eᴄ figures.

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `vae_models/`, `saved_sims/*.npy`
- **In this repo** (needed to reproduce the figures without a download): `saved_data/`, `params_I*.pkl` generation params

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

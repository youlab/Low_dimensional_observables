# Short / sparse time series

Truncated and sparsely-sampled versions of the gLV trajectories. Tests robustness of Eᴄ to limited temporal sampling.

## Model
`VAE_model.py`

## Run order
1. `train_truncate_sparse.py` / `train_sparse.py` — train VAEs on truncated / sparse series (`resubmit_truncate.sh`, `resubmit_sparse.sh`). Reads sims from `../VAE_gLV_simulation/`.
2. `Summary_short_sequences.ipynb` — FUV / Eᴄ figures.

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `vae_models/`
- **In this repo** (needed to reproduce the figures without a download): `saved_data/`

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

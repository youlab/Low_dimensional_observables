# Architecture benchmark

Compares VAE architecture variants against plain autoencoder and transformer baselines on the gLV data, to justify the model choice. Loads simulation data from `../../simulations/VAE_gLV_simulation/`.

## Model
`VAE_var2.py`, `VAE_var3.py`, `AutoEncoder_model.py`, `TransformerAE.py`

## Run order
1. `train_VAE_var2.py` / `train_VAE_var3.py` / `train_autoencoder.py` / `train_transformer.py` — train each architecture (`VAE_v2.sh`, `VAE_v3.sh`, `autoencoder_*.sh`, `transformer_*.sh`).
2. `Summary_VAE_var2.ipynb`, `Summary_VAE_var3.ipynb`, `Summary_AE.ipynb`, `Summary_transformer.ipynb` — per-model FUV.
3. `compare_VAEs.ipynb` / `benchmark.ipynb` — cross-architecture comparison figures.

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `AE_models/`, `VAE_v2/`, `VAE_v3/`, `transformer/`
- **In this repo** (needed to reproduce the figures without a download): `saved_data/` (reconstruction caches)

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

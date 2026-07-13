# Architecture benchmark

Compares the VAE against alternative compression models on the gLV data, to justify the
model choice: **VAE** (main), **VAE_v2**, **VAE_v3**, a plain **linear autoencoder**, a
**transformer** autoencoder, and **PCA** (a linear, non-learned baseline). Loads simulation
data from `../../simulations/VAE_gLV_simulation/`.

## Model
`VAE_var2.py`, `VAE_var3.py`, `AutoEncoder_model.py`, `TransformerAE.py` (PCA is computed
directly from the data — no trained model).

## Run order
1. `train_VAE_var2.py` / `train_VAE_var3.py` / `train_autoencoder.py` / `train_transformer.py`
   — train each architecture (`VAE_v2.sh`, `VAE_v3.sh`, `autoencoder_*.sh`, `transformer_*.sh`).
2. `Summary_VAE_var2.ipynb`, `Summary_VAE_var3.ipynb`, `Summary_AE.ipynb`,
   `Summary_transformer.ipynb` — per-model FUV.
3. `benchmark.ipynb` — cross-method comparison of FUV vs. embedding dimension. Three blocks,
   all at `model_index=1`: **bgLV fixed**, **bgLV random**, and **dgLV fixed**. Each block
   overlays all six methods per target-population size; PCA is the light-blue baseline.

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `AE_models/`, `VAE_v2/`, `VAE_v3/`, `transformer/`
- **In this repo** (needed to reproduce the figures without a download): `saved_data/`
  reconstruction caches for every method, including `saved_data/PCA/{bgLV,dgLV}/I*/` —
  the PCA FUV curves (single deterministic curve per target, no trials).

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

# Engineered target populations (fluorescent proteins)

VAE analysis of measured growth dynamics of engineered fluorescent-protein (FP) populations and community-assembly experiments.

## Model
`VAE_model.py`, `MLPVAE.py`

## Run order
1. `train_exp_target.py` / `train_growth_curve_inference.py` (+ `MLPVAE.py`) — train VAEs per FP target (`Train_EGFP_curve.sh`, `Train_mCherry_curve.sh`, `Train_mTagBFP2_curve.sh`, `Train_LSSmOrange_curve.sh`, `Train_comm_assembly.sh`).
2. `VAE_on_FP_populations.ipynb` / `VAE_on_single_target.ipynb` — Eᴄ / reconstruction figures. Use the in-repo `Nano_target_dynamics_{train,test}.npy` splits.

## Sequencing / community-composition analyses
Standalone scripts characterizing the measured community from the amplicon data in
`sequenced_data/` (independent of the VAE models — they run directly from the in-repo data):

- **`community_composition.py`** — normalizes `sequenced_data/sequence_composition.txt` to
  relative abundances, orders samples by Bray–Curtis hierarchical clustering, and plots the
  stacked strain composition with an inverse-Simpson diversity curve. Prints richness /
  diversity summary stats (289 samples).
  ```bash
  python community_composition.py
  ```
- **`infer_correlation.py`** — reconstructs absolute abundances
  (`sequence_composition.txt` × final OD from `sequenced_target_normal.npy` and
  `max_normalization.txt`), computes zero-filtered pairwise Pearson correlations, reports the
  frequency of significant pairs, and plots the correlation heatmap + per-target scatter fits.
  ```bash
  python infer_correlation.py
  ```
- **`confirm_composition_correspondence.py`** — checks which composition column belongs to
  which fluorescent protein, correlating each column against the named growth curves in
  `saved_data/growth_curves/<FP>/<FP>_ground_truth.npy`.
  ```bash
  python confirm_composition_correspondence.py
  ```
All three are verified to run under the project environment.

### Strain ordering
The 62 columns of `sequence_composition.txt` follow the strain-barcode order of the
amplicon pipeline, which for the four FP targets is **not** the manuscript K1–K4 order:
columns 0–3 are mCherry (K4), EGFP (K1), mTagBFP2 (K2), LSSmOrange (K3). Both analysis
scripts therefore apply `fp_col_order = [1, 2, 3, 0]` after loading, so k1–k4 in the
figures mean EGFP, mTagBFP2, LSSmOrange, mCherry. Strains 5–62 are unaffected. The
assignment is confirmed empirically by `confirm_composition_correspondence.py`.

## Data
- **On Zenodo** (restore with `bash ../../zenodo/zenodo_download.sh <DOI>`; see [`../../DATA.md`](../../DATA.md)): `vae_models/`, `mlp_models/`
- **In this repo** (needed to reproduce the figures without a download): `Nano_target_dynamics*.npy` (measured train/test splits), `sequenced_data/`

Notebook outputs are kept, so the published figures are visible without rerunning.
Figure-saving (`savefig`) is disabled for the release; plots still render inline.

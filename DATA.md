# Data & model archive

The trained models and heavy simulation arrays are **not** stored in this repository —
they are archived on Zenodo. The repository itself contains all source code, notebooks,
and the small derived results needed to reproduce the paper figures (FUV caches in
`*_embedding_FUV_CV/` and `saved_data/`, OLS-forecast outputs, `.mat` generation params,
and the experimental train/test splits).

**Zenodo DOI:** `10.5281/zenodo.XXXXXXX`  ← _fill in after publishing_

## Restore the artifacts

```bash
# from the repo root
bash zenodo/zenodo_download.sh 10.5281/zenodo.XXXXXXX            # everything
bash zenodo/zenodo_download.sh 10.5281/zenodo.XXXXXXX glv_models # a single tarball
```

Each tarball preserves repo-relative paths, so extraction restores every file to the
folder its notebook/script expects.

## Manifest

| Tarball | Restores into | Contents |
|---|---|---|
| `glv_models.tar.gz` | `simulations/VAE_gLV_simulation/{vae_models,mlp_models}/` | trained gLV VAE/MLP weights (`.pth`) |
| `glv_saved_sims.tar.gz` | `simulations/VAE_gLV_simulation/saved_sims/` | gLV simulation train/test arrays (`.npy`) |
| `glv_reconstruction.tar.gz` | `simulations/VAE_gLV_simulation/mlpvae_reconstruction/` | reconstruction arrays (`.npy`) |
| `benchmark_models.tar.gz` | `benchmarks/benchmark/{AE_models,VAE_v2,VAE_v3,transformer}/` | architecture-benchmark weights |
| `consumer_resource.tar.gz` | `simulations/consumer_resource/{vae_models,saved_sims}/` | MiCRM models + sims |
| `stochastic_logistic.tar.gz` | `simulations/stochastic_logistic/{vae_models,saved_sims}/` | models + sims |
| `network_effect.tar.gz` | `simulations/network_effect/{vae_models,saved_sims}/` | models + sims |
| `short_sequence.tar.gz` | `simulations/short_sequence/vae_models/` | truncated/sparse VAE weights |
| `lab_microbiome_models.tar.gz` | `experiments/lab_microbiome/{vae_models_CV,mlp_models}/` | CV VAE weights + collapse MLP weights |
| `lab_microbiome_datasets_CV.tar.gz` | `experiments/lab_microbiome/lab_microbiome_datasets_CV/` | CV datasets (1.3 GB) |
| `vaginal_microbiome.tar.gz` | `experiments/vaginal_microbiome/{vae_models_CV,vaginal_dataset_CV}/` | CV weights + datasets |
| `experimental_target.tar.gz` | `experiments/Keio_community/{vae_models,mlp_models}/` | FP-population model weights |
| `simulation_code_rawdata.tar.gz` | `simulation_code/network sparsity/` | raw generator example data (`random_init.txt`, 17.6 MB) |

Checksums for every tarball are in `SHA256SUMS` (uploaded alongside the archives and
verified automatically by `zenodo_download.sh`).

## What you do NOT need the archive for
The published FUV / embedding-dimension figures reproduce directly from the in-repo caches
(`*_embedding_FUV_CV/*.txt`, `saved_data/**/reconstruction_*.txt`). You only need the
Zenodo artifacts to **re-train** models or **re-compute** FUV/reconstructions from scratch.

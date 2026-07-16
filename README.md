# Low-dimensional observables of ecological dynamics

Code and analysis notebooks for learning **low-dimensional observables of ecological and
microbial community dynamics with variational autoencoders (VAEs)**. A VAE is trained to
compress high-dimensional community trajectories into a small latent embedding; the
minimum embedding dimension that reconstructs the dynamics (measured by the fraction of
unexplained variance, **FUV**) is the effective number of observables — the *Ec* — of the
system.

This repository accompanies the manuscript and reproduces every figure. It is organized as
self-contained folders, one per model system, grouped into simulations, architecture
benchmarks, and real-data experiments.

## Repository layout

```
simulations/          # synthetic dynamical systems the method is validated on
├── VAE_gLV_simulation/    generalized Lotka–Volterra (background + target communities)
├── consumer_resource/     MacArthur consumer-resource model (MiCRM)
├── stochastic_logistic/   stochastic logistic growth
├── network_effect/        interaction-network structure effects
└── short_sequence/        truncated / sparsely-sampled time series
benchmarks/
└── benchmark/         VAE architecture variants vs. autoencoder / transformer baselines
experiments/          # real measured communities
├── lab_microbiome/        soil & water lab communities (cross-validated)
├── vaginal_microbiome/    24-subject vaginal microbiome time series
└── Keio_community/        engineered fluorescent-protein populations (Keio strains)
simulation_code/      # MATLAB/Python generators for the synthetic dynamics (see its README)
```

Each system folder has its own `README.md` with the exact run order.

## Installation

```bash
conda env create -f environment.yml     # creates the `low_dim_observables` env
conda activate low_dim_observables
# or, with pip:  pip install -r requirements.txt   (Python 3.12)
```

Core dependencies: PyTorch 2.5, NumPy, SciPy, scikit-learn, pandas, matplotlib.

## Data & trained models

Source code and the **small derived results needed to reproduce the figures** live in this
repo (FUV caches, reconstruction summaries, OLS-forecast outputs, `.mat` generation params,
experimental train/test splits). The **trained model weights and heavy simulation arrays
are archived on Zenodo** — see [`DATA.md`](DATA.md).

- **Reproduce the published FUV / embedding-dimension figures:** no download needed — run
  the `Summary_*` notebooks, which read the in-repo caches.
- **Re-train models or recompute reconstructions from scratch:** restore the artifacts with
  `bash zenodo/zenodo_download.sh <DOI>` first.

## System → analysis map

| Folder | System | Produces |
|---|---|---|
| `simulations/VAE_gLV_simulation` | generalized Lotka–Volterra | Ec vs. community size; embedding-dimension scaling (main validation) |
| `simulations/consumer_resource` | consumer-resource (MiCRM) | Ec for resource-mediated dynamics |
| `simulations/stochastic_logistic` | stochastic logistic | Ec under intrinsic noise |
| `simulations/network_effect` | interaction networks | dependence of Ec on network structure |
| `simulations/short_sequence` | truncated/sparse series | robustness to short / sparse sampling |
| `benchmarks/benchmark` | architecture study | VAE variants vs. AE / transformer baselines |
| `experiments/lab_microbiome` | soil & water communities | Ec (CV), collapse prediction, community-shift & OLS forecasting |
| `experiments/vaginal_microbiome` | vaginal microbiome | Ec across 24 subjects |
| `experiments/Keio_community` | engineered FP populations (Keio strains) | Ec of measured growth dynamics |

> _Fill in the manuscript figure numbers for each row before release._

## Citing

If you use this code, please cite the manuscript and the Zenodo archive (DOI in
[`DATA.md`](DATA.md)). Licensed under the MIT License — see [`LICENSE`](LICENSE).

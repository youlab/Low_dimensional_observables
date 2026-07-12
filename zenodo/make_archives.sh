#!/usr/bin/env bash
# make_archives.sh — bundle the heavy artifacts (trained models + simulation arrays)
# into per-system tarballs for upload to Zenodo. Run from the repo root:
#
#     bash zenodo/make_archives.sh
#
# Output: zenodo/archives/*.tar.gz  and  zenodo/archives/SHA256SUMS
#
# What is EXCLUDED from the tarballs (these stay in the git repo):
#   *.mat, *.pkl        -> simulation generation params (kept in-repo for provenance)
#   *_traindoc.pth      -> tiny loss-doc arrays read by kept notebooks
# Each tarball preserves repo-relative paths, so `tar xzf <t> -C <repo-root>` restores
# every file to its original location (see zenodo/zenodo_download.sh).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
OUT="zenodo/archives"
mkdir -p "$OUT"

EXCLUDES=(--exclude='*.mat' --exclude='*.pkl' --exclude='*_traindoc.pth')

# tarball_name  <TAB>  space-separated source paths (repo-relative)
read -r -d '' GROUPS <<'EOF' || true
glv_models	simulations/VAE_gLV_simulation/vae_models simulations/VAE_gLV_simulation/mlp_models
glv_saved_sims	simulations/VAE_gLV_simulation/saved_sims
glv_reconstruction	simulations/VAE_gLV_simulation/mlpvae_reconstruction
benchmark_models	benchmarks/benchmark/AE_models benchmarks/benchmark/VAE_v2 benchmarks/benchmark/VAE_v3 benchmarks/benchmark/transformer
consumer_resource	simulations/consumer_resource/vae_models simulations/consumer_resource/saved_sims
stochastic_logistic	simulations/stochastic_logistic/vae_models simulations/stochastic_logistic/saved_sims
network_effect	simulations/network_effect/vae_models simulations/network_effect/saved_sims
short_sequence	simulations/short_sequence/vae_models
lab_microbiome_models	experiments/lab_microbiome/vae_models_CV experiments/lab_microbiome/mlp_models
lab_microbiome_datasets_CV	experiments/lab_microbiome/lab_microbiome_datasets_CV
vaginal_microbiome	experiments/vaginal_microbiome/vae_models_CV experiments/vaginal_microbiome/vaginal_dataset_CV
experimental_target	experiments/Keio_community/vae_models experiments/Keio_community/mlp_models
EOF

: > "$OUT/SHA256SUMS"
while IFS=$'\t' read -r name paths; do
  [ -z "$name" ] && continue
  # keep only paths that actually exist
  existing=()
  for p in $paths; do [ -e "$p" ] && existing+=("$p"); done
  if [ ${#existing[@]} -eq 0 ]; then echo "SKIP $name (no source paths present)"; continue; fi
  tarball="$OUT/${name}.tar.gz"
  echo "==> $tarball  <=  ${existing[*]}"
  tar "${EXCLUDES[@]}" -czf "$tarball" "${existing[@]}"
  ( cd "$OUT" && sha256sum "${name}.tar.gz" >> "SHA256SUMS" )
done <<< "$GROUPS"

# Raw generator data whose paths contain spaces (handled explicitly, outside the
# word-split GROUPS loop above).
RAWDATA=("simulation_code/network sparsity/random_init.txt")
present=(); for p in "${RAWDATA[@]}"; do [ -e "$p" ] && present+=("$p"); done
if [ ${#present[@]} -gt 0 ]; then
  echo "==> $OUT/simulation_code_rawdata.tar.gz  <=  ${present[*]}"
  tar -czf "$OUT/simulation_code_rawdata.tar.gz" "${present[@]}"
  ( cd "$OUT" && sha256sum "simulation_code_rawdata.tar.gz" >> "SHA256SUMS" )
fi

echo
echo "Done. Archives in $OUT/"
du -sh "$OUT"/*.tar.gz
echo "Checksums: $OUT/SHA256SUMS"

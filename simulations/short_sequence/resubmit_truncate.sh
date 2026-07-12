#!/bin/bash
# ============================================================
# resubmit_truncate.sh — retrain truncate models, save .pth
# Strategy: 1 node, 1 GPU, 4 trials in parallel per batch
# Submit 4 copies for 4-node parallelism:
#   sbatch resubmit_truncate.sh 1    ← node 1: model_idx=1
#   sbatch resubmit_truncate.sh 2    ← node 2: model_idx=2
#   sbatch resubmit_truncate.sh 3    ← node 3: model_idx=3
#   sbatch resubmit_truncate.sh 1 2  ← node 4: overflow/extras
# ============================================================
#SBATCH -p youlab-gpu
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=trunc_retrain
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kw367@duke.edu
#SBATCH --output=/hpc/group/youlab/wendyw/logs/slurm/trunc_retrain_%j.out

cd /hpc/group/youlab/wendyw/VAE_demo
source /opt/apps/rhel9/Anaconda3-2024.02/etc/profile.d/conda.sh
conda activate low_dim_observables

mkdir -p ./slurm_outputs/truncate
mkdir -p /hpc/group/youlab/wendyw/logs/slurm

# Which model indices to run — pass as args, default all
MODELS=${@:-"1 2 3"}
PARALLEL=4  # simultaneous trials on one GPU

echo "=============================="
echo "Node:    $SLURMD_NODENAME"
echo "GPU:     $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader)"
echo "Models:  $MODELS"
echo "Parallel: $PARALLEL"
echo "Start:   $(date)"
echo "=============================="

running=0

maybe_run() {
    local outfile=$1; shift
    # Skip if already successfully completed
    if [ -f "$outfile" ] && grep -q "Done:" "$outfile" 2>/dev/null; then
        return
    fi
    python -u "$@" > "$outfile" 2>&1 &
    running=$((running + 1))
    if (( running % PARALLEL == 0 )); then
        wait
        echo "  [$(date +%H:%M:%S)] batch complete (total launched: $running)"
    fi
}

for MODEL in $MODELS; do
  echo "--- model_idx=${MODEL} ---"
  for NT in 2 3 5 8 10; do
    for N_OBS in 5 10 20 50; do
      for TRIAL in 1 2 3 4 5; do
        OUT="./slurm_outputs/truncate_v2/trunc_M${MODEL}_NT${NT}_NOBS${N_OBS}_T${TRIAL}.out"
        maybe_run "$OUT" \
          train_truncate_sparse.py "${MODEL}" "100" "${NT}" "random" "${TRIAL}" "${N_OBS}"
      done
    done
  done
done

wait
echo "=============================="
echo "All done: $(date)"
echo "Truncate .pth files saved:"
ls ./vae_models/bgLV/I*/truncate/truncate_*.pth 2>/dev/null | wc -l
echo "=============================="
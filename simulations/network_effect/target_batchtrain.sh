#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --ntasks=16
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=3-00:00:00
#SBATCH --exclude=dcc-youlab-gpu-ferc-s-aa25-19

set -euo pipefail

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

mkdir -p ./slurm_outputs

MAX_JOBS=16

MODEL=5
for TRIAL in 4 5; do
  echo "=== TRIAL ${TRIAL} ==="
  for TGID in {0..100}; do
    echo "--- TGID ${TGID} ---"

    srun --nodes=1 \
        --ntasks=1 \
        --gpus-per-task=1 \
        --cpus-per-task=1 \
        --exact \
        --output="./slurm_outputs/M${MODEL}_TG${TGID}_T${TRIAL}.out" \
        --error="./slurm_outputs/M${MODEL}_TG${TGID}_T${TRIAL}.err" \
        python -u train_indiv_target.py "${TGID}" "${MODEL}" "${TRIAL}" &

    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
      sleep 5
    done
  done
done

wait
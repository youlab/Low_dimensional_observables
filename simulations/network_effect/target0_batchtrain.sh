#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --ntasks=3
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

MODEL=5
TGID=0

for TRIAL in 1 2 3; do
  echo "=== TRIAL ${TRIAL} ==="
  echo "--- TGID ${TGID} ---"

  srun --nodes=1 \
      --ntasks=1 \
      --gpus-per-task=1 \
      --cpus-per-task=1 \
      --exact \
      --output="./slurm_outputs/M${MODEL}_TG${TGID}_T${TRIAL}.out" \
      --error="./slurm_outputs/M${MODEL}_TG${TGID}_T${TRIAL}.err" \
      python -u train_indiv_target.py "${TGID}" "${MODEL}" "${TRIAL}" &


done

wait
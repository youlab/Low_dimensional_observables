#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --ntasks=9
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2-00:00:00
#SBATCH --exclude=dcc-youlab-gpu-ferc-s-aa25-19

set -euo pipefail

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# Create output folder if it doesn't exist
mkdir -p ./slurm_outputs

for MODEL in 2 3 4 5; do
  echo "=== MODEL ${MODEL} ==="
  for TRIAL in 1 2 3 4 5; do
    echo "=== TRIAL ${TRIAL} ==="
      for NT in 2 3 4 5 6 7 8 9 10; do
      echo "--- NT ${NT} ---"
      srun --nodes=1 \
          --ntasks=1 \
          --gpus-per-task=1 \
          --cpus-per-task=1 \
          --exact \
          --output="./slurm_outputs/network_M${MODEL}_NT${NT}_T${TRIAL}.out" \
          --error="./slurm_outputs/network_M${MODEL}_NT${NT}_T${TRIAL}.err" \
          python -u train_network_gLVs.py "${NT}" "${MODEL}" "${TRIAL}" &
      done
  done
  wait
done
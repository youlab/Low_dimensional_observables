#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --job-name=vag_vae_cv
#SBATCH --nodes=8                   # request 8 nodes
#SBATCH --ntasks=8                  # 8 tasks (1 per node)
#SBATCH --gres=gpu:1                # 1 GPU per node
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclusive

# Cross-validated VAE training for one n_target.
# Usage: sbatch train_vag_CV.sh <N_TARGET>
#   e.g. sbatch train_vag_CV.sh 1
# Sweeps the 24 held-out subjects x 3 training trials; each run covers all
# embedding dimensions in Train_vag_microbiome_CV.py.

if [ "$#" -ne 1 ]; then
  echo "Usage: sbatch $0 <N_TARGET>   # N_TARGET in 1..5" >&2
  exit 1
fi

NT=$1

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# Create output folders if they don't exist
mkdir -p ./slurm_outputs ./vae_models_CV

# One srun per held-out subject, launched in waves of 8 (one per node)
for TRIAL in 1 2 3; do
  for SUB in $(seq 1 24); do
    srun --ntasks=1 -N1 --exclusive \
      --output=./slurm_outputs/NT${NT}_sub${SUB}_trial${TRIAL}.out \
      python Train_vag_microbiome_CV.py $NT $SUB $TRIAL &
    if (( SUB % 8 == 0 )); then
      wait  # keep at most 8 concurrent tasks
    fi
  done
  wait
done

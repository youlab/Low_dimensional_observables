#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --job-name=lab_vae_cv
#SBATCH --nodes=8                   # request 8 nodes
#SBATCH --ntasks=8                  # 8 tasks (1 per node)
#SBATCH --gres=gpu:1                # 1 GPU per node
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclusive

# Cross-validated VAE training for one (community, n_target) pair.
# Usage: sbatch train_lab_CV.sh <SAMPLE> <N_TARGET>
#   e.g. sbatch train_lab_CV.sh Soil-A 1
# Sweeps the 8 held-out replicates x 3 training trials; each run covers all
# embedding dimensions in Train_lab_microbiome_CV.py.

if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch $0 <SAMPLE> <N_TARGET>" >&2
  echo "  SAMPLE in {Soil-A, Soil-B, Soil-C, Water-A, Water-B}; N_TARGET in 1..5" >&2
  exit 1
fi

SAMPLE=$1
NT=$2

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# Create output folders if they don't exist
mkdir -p ./slurm_outputs ./vae_models_CV

# One srun per held-out replicate; 8 run concurrently, one trial at a time
for TRIAL in 1 2 3; do
  for REP in 1 2 3 4 5 6 7 8; do
    srun --ntasks=1 -N1 --exclusive \
      --output=./slurm_outputs/${SAMPLE}_NT${NT}_rep${REP}_trial${TRIAL}.out \
      python Train_lab_microbiome_CV.py $SAMPLE $NT $REP $TRIAL &
  done
  wait  # finish all 8 folds of this trial before starting the next
done

#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --nodes=10                   # request 5 nodes
#SBATCH --ntasks=10                  # 5 tasks (1 per node)
#SBATCH --gres=gpu:1                # 1 GPU per node
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclusive

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

FP="mTagBFP2"

# Create output folder if it doesn't exist
mkdir -p ./slurm_outputs
mkdir -p ./mlp_models/${FP}

# Use srun to launch each job on a different task/node
for TRIAL in 1 2 3 4 5 6 7 8 9 10; do
  srun --gres=gpu:1 --ntasks=1 -N1 --exclusive \
    --output=./slurm_outputs/${FP}_T${TRIAL}.out \
    python train_growth_curve_inference.py $FP $TRIAL &
done
wait  # Wait for all 5 jobs before launching next batch

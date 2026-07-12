#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --nodes=5                   # request 5 nodes
#SBATCH --ntasks=5                  # 5 tasks (1 per node)
#SBATCH --gres=gpu:1                # 1 GPU per node
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclusive

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# Create output folder if it doesn't exist
mkdir -p ./slurm_outputs

# Use srun to launch each job on a different task/node
for TRIAL in 1 2 3 4 5; do
  for NT in 1 2 3 4 5; do
    srun --ntasks=1 -N1 --exclusive \
      --output=./slurm_outputs/train_NT${NT}_T${TRIAL}.out \
      python train_exp_target.py $NT $TRIAL &
  done
  wait  # Wait for all 5 jobs before launching next batch
done

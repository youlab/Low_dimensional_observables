#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --nodes=5                   # request 5 nodes
#SBATCH --ntasks=5                  # 5 tasks (1 per node)
#SBATCH --gres=gpu:1                # 1 GPU per node
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=./slurm_outputs/retrain_full.out
#SBATCH --exclusive

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# Use srun to launch each job on a different task/node
srun --ntasks=1 --exclusive -N1 -n1 python Retrain_lab_microbiome.py "Soil-A" &
srun --ntasks=1 --exclusive -N1 -n1 python Retrain_lab_microbiome.py "Soil-B" &
srun --ntasks=1 --exclusive -N1 -n1 python Retrain_lab_microbiome.py "Soil-C" &
srun --ntasks=1 --exclusive -N1 -n1 python Retrain_lab_microbiome.py "Water-A" &
srun --ntasks=1 --exclusive -N1 -n1 python Retrain_lab_microbiome.py "Water-B" &
wait

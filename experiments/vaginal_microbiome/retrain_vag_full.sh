#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --job-name=retrain_Ec
#SBATCH --output=./slurm_outputs/retrain_all.out
#SBATCH --exclusive

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

mkdir -p ./slurm_outputs

srun python Retrain_vag_microbiome_Ec.py 5 &
srun python Retrain_vag_microbiome_Ec.py 8 &
srun python Retrain_vag_microbiome_Ec.py 10 &
wait

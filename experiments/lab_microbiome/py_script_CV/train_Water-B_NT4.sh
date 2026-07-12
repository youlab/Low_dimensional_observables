#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --nodes=8                   # request 8 nodes
#SBATCH --ntasks=8                  # 8 tasks (1 per node)
#SBATCH --gres=gpu:1                # 1 GPU per node
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=./slurm_outputs/Water-B_NT4.out
#SBATCH --exclusive

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

SAMPLE="Water-B"
NT=4

# Use srun to launch each job on a different task/node
srun --ntasks=1 --exclusive -N1 -n1 python Train_lab_microbiome_CV.py $SAMPLE $NT 1 &
srun --ntasks=1 --exclusive -N1 -n1 python Train_lab_microbiome_CV.py $SAMPLE $NT 2 &
srun --ntasks=1 --exclusive -N1 -n1 python Train_lab_microbiome_CV.py $SAMPLE $NT 3 &
srun --ntasks=1 --exclusive -N1 -n1 python Train_lab_microbiome_CV.py $SAMPLE $NT 4 &
srun --ntasks=1 --exclusive -N1 -n1 python Train_lab_microbiome_CV.py $SAMPLE $NT 5 &
srun --ntasks=1 --exclusive -N1 -n1 python Train_lab_microbiome_CV.py $SAMPLE $NT 6 &
srun --ntasks=1 --exclusive -N1 -n1 python Train_lab_microbiome_CV.py $SAMPLE $NT 7 &
srun --ntasks=1 --exclusive -N1 -n1 python Train_lab_microbiome_CV.py $SAMPLE $NT 8 &
wait

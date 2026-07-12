#!/bin/bash
#SBATCH -p youlab-gpu
#SBATCH --job-name=collapse_mlp
#SBATCH --nodes=5                   # request 5 nodes
#SBATCH --ntasks=5                  # 5 tasks (1 per node)
#SBATCH --gres=gpu:1                # 1 GPU per node
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclusive
#SBATCH --exclude=dcc-youlab-gpu-28

# Load environment
source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# Create output folder if it doesn't exist
mkdir -p ./slurm_outputs

THRESHOLD=0.2

# Use srun to launch each job on a different task/node
for SAMPLE in "Soil-A" "Soil-B" "Soil-C" "Water-A" "Water-B"; do
  for NASV in 1 2 3 4 5 6 7 8 9 10; do
    for TRIAL in 1 2 3 4 5; do
      srun --ntasks=1 -N1 --exclusive \
        --output=./slurm_outputs/collapse_${SAMPLE}_${NASV}_${THRESHOLD}_T${TRIAL}.out \
        python MLP_collapse_prediction.py $SAMPLE $NASV $THRESHOLD $TRIAL &
    done
    wait  # Wait for all 5 jobs before launching next batch
  done
done

#!/bin/bash
#SBATCH --job-name=Nf9
#SBATCH --array=1-10000
#SBATCH -p youlab
#SBATCH --mem=4G

module load Matlab/R2020a
matlab -nodesktop -nodisplay -singleCompThread -r "rank=$SLURM_ARRAY_TASK_ID; sgLV_v04a; quit"

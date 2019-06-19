#!/bin/bash
#SBATCH -c 4
#SBATCH --mem 32G
#SBATCH --mail-user=a.solovyev39@gmail.com
#SBATCH -o slurm-%J.out
#SBATCH -e slurm-%J.err
#SBATCH --array 1-40%10
Rscript --vanilla create_hm_matrix.R "A549"
echo "binned signal calculated"


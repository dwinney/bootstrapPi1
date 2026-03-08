#! /usr/bin/bash

#SBATCH --job-name=pi1_BS
#SBATCH --output=output.txt
#SBATCH --partition=icn
#SBATCH --ntasks=1
#SBATCH --nodes=1

srun ./do_fit
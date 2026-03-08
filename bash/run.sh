#! /usr/bin/bash

#SBATCH --job-name=pi1_BS
#SBATCH --output=output.txt
#SBATCH --partition=icn
#SBATCH --ntasks=1
#SBATCH --nodes=1

module purge
module load lamod/boost/1.80
module load lamod/ROOT/v6-28-04
srun ./do_fit
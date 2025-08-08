#!/bin/bash -e

#SBATCH --job-name=05-PD-view-cov
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=00:05:00
#SBATCH --mem=28G
#SBATCH --ntasks=1
#SBATCH --profile=task
#SBATCH --account=ga03048
#SBATCH --cpus-per-task=2

# Takes two arguments: 1) PRI/ALT and 2) R1/R2

bash ./05-PD-view-cov.sh PRI R1

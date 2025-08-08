#!/bin/bash -e

#SBATCH --job-name=04-purge-dups
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=01:20:00
#SBATCH --mem=38G
#SBATCH --ntasks=1
#SBATCH --profile=task
#SBATCH --account=ga03048
#SBATCH --cpus-per-task=2

# Takes two arguments: 1) PRI/ALT and 2) R1/R2

bash ./04-PD-purge.sh PRI R1

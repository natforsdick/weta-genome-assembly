#!/bin/bash -e

#SBATCH --job-name=make-paf
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=22G
#SBATCH --ntasks=1
#SBATCH --profile=task 
#SBATCH --account=ga03048
#SBATCH --cpus-per-task=62

# Takes one param: PRI/ALT

bash ./01-PD-make-paf.sh PRI

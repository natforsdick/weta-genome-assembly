#!/bin/bash -e

#SBATCH --job-name=minimap-aln
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=120G #Used 26 for whole genomes
#SBATCH --profile=task
#SBATCH --account=ga03048
#SBATCH --cpus-per-task=8 # Used 32 for whole genomes

bash ./01-align-genomes.sh

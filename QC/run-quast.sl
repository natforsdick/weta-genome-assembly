#!/bin/bash -e
#SBATCH -J quast
#SBATCH -A ga03048
#SBATCH --cpus-per-task 8
#SBATCH --time 00:25:00
#SBATCH --mem=16G
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

bash ./quast.sh

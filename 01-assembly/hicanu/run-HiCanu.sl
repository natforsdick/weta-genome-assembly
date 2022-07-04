#!/bin/bash -e

#SBATCH --account ga03048
#SBATCH --job-name hicanu
#SBATCH --mem 100M
#SBATCH --time 00:01:00
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err

bash ./run-HiCanu.sh

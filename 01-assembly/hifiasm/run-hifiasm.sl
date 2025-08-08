#!/bin/bash -e

#SBATCH --account ga03048
#SBATCH --job-name hifiasm
#SBATCH --cpus-per-task=64 # When starting this job, set to 24 CPU, then increase to match mem when OOM
#SBATCH --mem 250G # When initially starting this job, set to 22 GB. Then when OOM, restart with 50 GB
#SBATCH --time 03:00:00 # When starting this job, set to 04:00:00, then change to 00:30:00 for final mem-intensive step
#SBATCH --partition=bigmem
#SBATCH --output hifiasm.%j.out
#SBATCH --error hifiasm.%j.err
#SBATCH --profile=task

bash run-hifiasm.sh

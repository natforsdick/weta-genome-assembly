#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J nanoqc # job name (shows up in the queue)
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --array 1-10
#SBATCH --mem=100M
#SBATCH --partition=large
#SBATCH --time=00:10:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%A.out #
#SBATCH --error %x.%A.err #

# This script takes one variable $1 : /path/to/input/data/

bash ./nanoQC.sh /path/to/input/data/

#!/bin/bash -e

#SBATCH -A ga03048
#SBATCH -J BUSCO # job name (shows up in the queue)
#SBATCH -c 48
#SBATCH --mem=70GB
#SBATCH --partition=large
#SBATCH --time=05:00:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out #
#SBATCH --error %x.%j.err #

bash ./BUSCO3.sh

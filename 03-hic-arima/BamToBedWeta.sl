#!/bin/bash
#SBATCH --account=ga03048
#SBATCH --job-name=BamToBedWeta # job name (shows up in the queue)
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=01:15:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

bash ./BamToBedWeta.sh

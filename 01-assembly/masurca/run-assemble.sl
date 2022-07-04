#!/bin/bash -e

#SBATCH --account ga03048
#SBATCH --job-name masurca-weta
#SBATCH --cpus-per-task=48
#SBATCH --mem 140G 
#SBATCH --time 12:00:00
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

export TMPDIR=/nesi/nobackup/ga03048/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

bash ./run-assemble2.sh

#!/bin/bash -e

#SBATCH --account ga03048
#SBATCH --job-name masurca-weta
#SBATCH --cpus-per-task=12
#SBATCH --mem 50G 
#SBATCH --time 00:20:00
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

#############
# run-masurca.sl
# Nat Forsdick, 2021-11-01
# Masurca hybrid genome assembly
#############

############
# MODULES
ml purge
ml MaSuRCA/4.0.9-gimkl-2020a
ml minimap2/2.20-GCC-9.2.0
############

export TMPDIR=/nesi/nobackup/ga03186/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

echo Beginnning MaSuRCA at
date

cd /nesi/nobackup/ga03048/assemblies/masurca4.0.9/

bash ./assemble.sh > 2022-04-22-assemble.out

echo Finishing MaSuRCA at
date


#!/bin/bash -e

#SBATCH --account ga03048
#SBATCH --job-name masurca-weta
#SBATCH --cpus-per-task=64
#SBATCH --mem 96G #16G #32G
#SBATCH --time 12:00:00
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
ml MaSuRCA/4.0.5-gimkl-2020a
ml minimap2/2.20-GCC-9.2.0
############

echo Beginnning MaSuRCA at
date

cd /nesi/nobackup/ga03048/assemblies/masurca/

bash /nesi/nobackup/ga03048/assemblies/masurca/assemble.sh > 01-assemble.out

echo Finishing MaSuRCA at
date


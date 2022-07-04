#!/bin/bash -e

#############
# Assembly with HiCanu
# Nat Forsdick, 2022-03-09

#############
# MODULES
module purge
module load Canu/2.1.1-GCC-9.2.0
ml list
#############

date

# -d = assembly directory
canu -p weta-canu -d /nesi/nobackup/ga03048/assemblies/weta-canu genomeSize=6g \
batMemory=512 batThreads=32 minThreads=8 minMemory=16 \
useGrid=true gridOptions="--account=ga03048 --job=hicanu --time=04:00:00 --output=%x.%j.out --error=%x.%j.err" \
-pacbio /nesi/project/ga03048/data/pacbio/hifi/*.fastq.gz

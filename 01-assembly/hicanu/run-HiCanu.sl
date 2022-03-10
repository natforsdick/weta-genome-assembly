#!/bin/bash -e

#SBATCH --account ga03048
#SBATCH --job-name hicanu
#SBATCH --mem 100M
#SBATCH --time 00:01:00
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err


module purge
module load Canu/2.1.1-GCC-9.2.0
ml list

date

# -d = assembly directory
canu -p weta-canu -d /nesi/nobackup/ga03048/assemblies/weta-canu genomeSize=6g \
batMemory=512 batThreads=32 minThreads=8 minMemory=16 \
useGrid=true gridOptions="--account=ga03048 --job=hicanu --time=04:00:00 --output=%x.%j.out --error=%x.%j.err" \
-pacbio /nesi/project/ga03048/data/pacbio/hifi/*.fastq.gz

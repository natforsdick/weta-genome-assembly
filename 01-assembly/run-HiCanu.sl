#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J hicanu 
#SBATCH -n 1
#SBATCH --mem=512M
#SBATCH --partition=large
#SBATCH --time=00:05:00 
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err 

module purge
module load Canu/2.1.1-GCC-9.2.0
ml list

date

# -d = assembly directory
canu -p weta-canu -d /nesi/nobackup/ga03048/assemblies/weta-canu genomeSize=6g \
useGrid=true gridOptions="--account=ga03048 --job=hicanu --time=01:00:00 --output=%x.%j.out --error=%x.%j.err --profile=task" \
-pacbio /nesi/project/ga03048/data/pacbio/hifi/*.fastq.gz

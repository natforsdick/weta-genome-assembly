#!/bin/bash -e 
#SBATCH -A ga03048
#SBATCH -J meryl-merge
#SBATCH -c 28
#SBATCH --mem=6G
#SBATCH --time=00:10:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

##########
# PARAMS #
##########
genome=weta
outdir=/nesi/nobackup/ga03048/assemblies/hifiasm/05-merqury/

##########

cd $outdir

echo "merging"
date
# 2. Merge
meryl union-sum threads=22 memory=$SLURM_MEM_PER_NODE output ${genome}.meryl ${outdir}*.meryl

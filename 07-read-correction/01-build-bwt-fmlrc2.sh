#!/bin/bash -e

#SBATCH -A ga03048
#SBATCH -J fmlrc2
#SBATCH --time 1-00:00:00 # step 1 used 20 hrs
#SBATCH --mem 80G # step 1: 40 GB
#SBATCH --cpus-per-task 24 # step 1: ~2 CPU
#SBATCH	--error=%x.%j.err
#SBATCH	--out=%x.%j.out
#SBATCH	--profile=task

# 01-build-bwt-fmlrc2.sh
# Nat Forsdick, 2022-04-06
# Step 1 in fmlrc2 correction of raw PacBio reads

###########
# PARAMS  #
datadir=/nesi/project/ga03048/weta/fastq/
illuminadir=/nesi/nobackup/ga03048/correction/trimmomatic/
outdir=/nesi/nobackup/ga03048/correction/fmlrc2/
export TMPDIR=/nesi/nobackup/ga03048/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR
###########

###########
# MODULES #
module purge
module load Anaconda3
source activate ropebwt2
module load rust-fmlrc/0.1.5-GCCcore-9.2.0
###########

echo "Beginning building BWT at"
date

cd ${outdir}

# step 1: make BWT
gunzip -c ${illuminadir}H07456-L1_S1_L001_R1_001*P.fastq.gz | awk "NR % 4 == 2" | sort -T $TMPDIR | tr NT TN |\
	ropebwt2 -LR | tr NT TN | fmlrc2-convert weta_msbwt.npy

echo "Completed BWT at"
date

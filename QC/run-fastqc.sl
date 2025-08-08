#!/bin/bash -e
#SBATCH -A ga03048 # CHANGE!!
#SBATCH -J fastqc # job name (shows up in the queue)
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH --partition=large
#SBATCH --time=00:40:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output=fastqc.%j.out
#SBATCH --error=fastqc.%j.err

###################
# FastQC 
# Nat Forsdick, 2020-12-10
###################

# This script takes two arguments, $1 : the path to the input directory, \
# and $2 : the prefix of the fastq file not including the lane number.

bash fastqc.sh /nesi/nobackup/ga03048/data/illumina/ H07456-L1_S1_L00*

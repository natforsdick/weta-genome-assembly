#!/bin/bash -e 
#SBATCH -A ga03048 
#SBATCH -J fastqc # job name (shows up in the queue) 
#SBATCH -N 1
#SBATCH -c 3
#SBATCH -n 1
#SBATCH --mem=900M 
#SBATCH --array=1-8
#SBATCH --partition=large 
#SBATCH --time=02:30:00 #Walltime (HH:MM:SS) 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz 
#SBATCH --output=Array_fastqc.%A.out
#SBATCH --error=Array_fastqc.%A.err     

################### 
# FastQC (as an array)
# Nat Forsdick, 2020-12-10
################### 

# Before running, be sure to adjust the array numbers to match the L00* of inputs, \
# and the runtime based on file sizes. F+R fastqs of ~60G take ~2.5 hrs to run.

# This script takes two arguments, $1 : the path to the input directory, \ 
# and $2 : the prefix of the fastq file not including the lane number.
# e.g. to execute:
# sbatch run_array_fastqc.sl /nesi/nobackup/ga03048/data/illumina/ H07456-L1_S1_L00

#### MODULES ####
module purge
module load FastQC/0.11.9 MultiQC/1.9-gimkl-2020a-Python-3.8.2
#################

#### ENVIRONMENT #### 
# Example sampname: H07456-L1_S1_L002_R1_001.fastq.gz  
# We pass the L00* in as the array numbers

IN_DIR=$1
#IN_PREFIX=$2

#fq1=_R1_001.fastq.gz
#fq2=_R2_001.fastq.gz
#####################

cd ${IN_DIR}

#echo "Running FastQC for ${SLURM_ARRAY_TASK_ID}"
#fastqc -t 2 ${IN_DIR}${IN_PREFIX}${SLURM_ARRAY_TASK_ID}${fq1} 
#fastqc -t 2 ${IN_DIR}${IN_PREFIX}${SLURM_ARRAY_TASK_ID}${fq2}
#echo "Completed FastQC for ${SLURM_ARRAY_TASK_ID}"

echo "Starting FastQC Job"
file=$(ls ./*.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
fastqc ${file}
echo "Finishing FastQC Job"

# Once finished, you can run MultiQC
#multiqc ${IN_DIR}

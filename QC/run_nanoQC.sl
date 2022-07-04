#!/bin/bash -e 
#SBATCH -A ga03048 
#SBATCH -J nanoqc # job name (shows up in the queue) 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 2
#SBATCH --array 1-10
#SBATCH --mem=100M 
#SBATCH --partition=large 
#SBATCH --time=00:10:00 #Walltime (HH:MM:SS) 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz 
#SBATCH --output %x.%A.out #  
#SBATCH --error %x.%A.err #     

################### 
# NanoQC (for PoreC/MinION data)
# Nat Forsdick  - 2020-12-03
################### 

# Running this to explore the read quality of the Hi-C data for weta and kaki.
# To be interpreted alongside Hi-C_QC reports.
# This script takes one variable $1 : /path/to/input/data/
# This only works if the directory only contains fastq files.
# To run:
# sbatch run_nanoQC.sl /path/to/fastqs/

#### MODULES #### 
module purge
module load Python/3.7.3-gimkl-2018b 
################### 

#### ENVIRONMENT ####
DATADIR=$1

cd ${DATADIR}
file=$(ls ${DATADIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)
#####################

nanoQC -o ${DATADIR}${file}_QC ${DATADIR}${file} 


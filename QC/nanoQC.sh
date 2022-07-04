#!/bin/bash -e 

############
# NanoQC (for PoreC/MinION data)
# Nat Forsdick  - 2020-12-03 

# Running this to explore the read quality of the Pore-C data for weta.
# To be interpreted alongside Hi-C_QC reports.
# This script takes one variable $1 : /path/to/input/data/
# This only works if the directory only contains fastq files.

############
# MODULES 
module purge
module load Python/3.7.3-gimkl-2018b 
############ 

# ENVIRONMENT
DATADIR=$1

cd ${DATADIR}
file=$(ls ${DATADIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)
###########

nanoQC -o ${DATADIR}${file}_QC ${DATADIR}${file} 


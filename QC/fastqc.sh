#!/bin/bash -e 

################### 
# FastQC
# Nat Forsdick, 2020-12-10
################### 

# This script takes two arguments, $1 : the path to the input directory, \ 
# and $2 : the prefix of the fastq file not including the lane number.
# e.g. to execute:
# bash run_fastqc.sh /nesi/nobackup/ga03048/data/illumina/ H07456-L1_S1_L00

#### MODULES ####
module purge
module load FastQC/0.11.9  
## MultiQC/1.9-gimkl-2020a-Python-3.8.2
#################

#### ENVIRONMENT #### 
# Example sampname: H07456-L1_S1_L002_R1_001.fastq.gz  

IN_DIR=$1
IN_PREFIX=$2

#fq1=_R1_001.fastq
#fq2=_R2_001.fastq
fq1=_1P.fq.gz
fq2=_2P.fq.gz
#####################

cd ${IN_DIR}

fastqc -t 4 ${IN_DIR}${IN_PREFIX}${fq1} 
fastqc -t 4 ${IN_DIR}${IN_PREFIX}${fq2}

# Once finished, you can run MultiQC
#multiqc ${IN_DIR}

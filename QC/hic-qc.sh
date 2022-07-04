#!/bin/bash -e 

################### 
# Run HiC_QC
# Nat Forsdick, 2020-11-27
################### 

# This script takes two arguments, $1 : the path to the input directory, / 
# and $2 : the prefix of the input file. 
# e.g. to execute: 
# sbatch run_hic_qc.sl /nesi/nobackup/ga03048/results/hic-qc-kaki/ Kaki_HiC_mapped 

################### 
# Need this version of miniconda
module purge
module load Miniconda3/4.7.10
################### 

#### ENVIRONMENT #### 
source /opt/nesi/CS400_centos7_bdw/Miniconda3/4.4.10/etc/profile.d/conda.sh
conda activate hic_qc

hic_qc=/nesi/project/ga03048/modules/hic_qc/hic_qc.py
IN_DIR=$1
IN_BAM=$2
#####################


cd $IN_DIR
echo "running hic_qc for ${IN_BAM}"
python ${hic_qc} -b ${IN_BAM}.bam -r -o ${IN_BAM}.hicqc 
echo "finished running hic_qc for ${IN_BAM}"

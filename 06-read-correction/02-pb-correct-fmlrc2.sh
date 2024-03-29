#!/bin/bash -e

# 02-pb-correct-fmlrc2.sh
# Nat Forsdick, 2022-04-07
# Step 2 for correcting raw pacbio CLR data
# Run 01-build-bwt-fmlrc2.sh first

###########
# PARAMS  #
datadir=/nesi/project/ga03048/weta/fastq/ # raw PacBio reads
illuminadir=/nesi/nobackup/ga03048/correction/trimmomatic/
outdir=/nesi/nobackup/ga03048/correction/fmlrc2/
###########

###########
# MODULES #
module purge
module load Anaconda3
source activate ropebwt2
module load rust-fmlrc/0.1.7-GCCcore-11.2.0
###########

# Avoid SLURM OOM issues
export TMPDIR=/nesi/nobackup/ga03186/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

cd ${datadir}

echo "Beginning fmlrc2 correction at"
date

#step 2: run fmlrc2 correction - $datadir holds raw PacBio reads

file=$(ls *.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)
filename=$(basename "$file")
filename=${filename%.*}

echo $filename

fmlrc2 -t $SLURM_CPUS_PER_TASK ${outdir}weta_msbwt.npy ${file} ${outdir}${filename}-fmlrc2-corr-reads.fa

echo "Completed fmlrc2 correction for ${filename} at"
date

#!/bin/bash -e
#SBATCH --account ga03048
#SBATCH --job-name 01-make-meryl
#SBATCH --time 00:30:00 # could need a couple of hours per fastq
#SBATCH --mem=130G # will prob need at least 24
#SBATCH --cpus-per-task=46
#SBATCH --output %j.%A.%a.out
#SBATCH --error %j.%A.%a.err
#SBATCH --profile=task
#SBATCH --array=1-2

##########
# PARAMS #
########## 
# to get best kmers, run:
# sh $MERQURY/best_k.sh [genome size in bp]
k=21 
genome=weta # output prefix

indir=/nesi/nobackup/ga03048/correction/trimmomatic/
outdir=/nesi/nobackup/ga03048/assemblies/hifiasm/05-merqury/
#########

cd $indir

# call samplist from file, and pass to array
SAMPLE_LIST=($(<input-fastq-list3.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]+1}

# test for single samp uses for loop:
#for file in ${indir}H07456-L1_S1_L001_R1_001.fastq_trimmed_1P.fastq.gz
#do

filename=$(basename "$SAMPLE") 
filename=${filename%.*} 
filename=${filename%.*}
echo $filename
date

meryl threads=32 memory=$SLURM_MEM_PER_NODE k=$k count output ${outdir}${filename}.meryl ${indir}${filename}.fastq.gz
date


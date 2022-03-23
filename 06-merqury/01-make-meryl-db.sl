#!/bin/bash -e
#SBATCH --account ga03048
#SBATCH --job-name 01-make-meryl
#SBATCH --time 00:02:00 # could need a couple of hours per fastq
#SBATCH --mem=200M #140G # will prob need at least 24
#SBATCH --cpus-per-task=2#46
#SBATCH --output %A.%a.out
#SBATCH --error %A.%a.err
#SBATCH --profile=task
#SBATCH --array=1-6

k=21
genome=weta
indir=/nesi/nobackup/ga03048/correction/trimmomatic/
outdir=/nesi/nobackup/ga03048/assemblies/hifiasm/05-merqury/


cd $indir

SAMPLE_LIST=($(<input-fastq-list.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

# loop for whole data set
#for file in `ls *P.fastq.gz`
#for file in ${indir}H07456-L1_S1_L001_R1_001.fastq_trimmed_1P.fastq.gz
#do

filename=$(basename "$SAMPLE") 
filename=${filename%.*} 
filename=${filename%.*}
echo $filename
date

echo "meryl threads=32 memory=$SLURM_MEM_PER_NODE k=$k count output ${outdir}${filename}.meryl ${indir}${filename}.fastq.gz"


#echo "merging"
#date
# 2. Merge
#meryl union-sum threads=32 memory=$SLURM_MEM_PER_NODE output ${genome}.meryl ${outdir}*.meryl

#!/bin/bash -e
#SBATCH --account=ga03186
#SBATCH --job-name=fastp 
#SBATCH --cpus-per-task=12 
#SBATCH --mem=16G
#SBATCH --time=6:00:00 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

##########
# PARAMS
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/
ASSEMBLY=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/02-purge-dups/02-weta-hic-hifiasm-p_ctg-purged.fa
HIC_DIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/
OUTDIR=2023-OmniC-out
READ1=_R1_001
READ2=_R2_001
fq=.fastq.gz
############

ml purge && module load fastp/0.23.2-GCC-11.3.0

cd ${HIC_DIR}$OUTDIR
### Clean HiC Reads with fastp.###
for file in ${HIC_DIR}*${READ1}${fq}
do

echo processing $file
f=${file%${READ1}${fq}}
y=${f##*/}

fastp \
-i ${file} \
-o ${y}${READ1}_clean${fq} \
-I ${HIC_DIR}${y}${READ2}${fq} \
-O ${y}${READ2}_clean${fq} \
--trim_front1 15 \
--trim_front2 15 \
--qualified_quality_phred 20 \
--length_required 50 \
--thread 12

done

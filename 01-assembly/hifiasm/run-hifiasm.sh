#!/bin/bash -e

##############
# HIFIASM WITH HI-C
# 2021-08-16
# Nat Forsdick
##############

##############
# MODULES
module purge
module load hifiasm/0.15.5-GCC-9.2.0
##############

##############
INDIR=/nesi/project/ga03048/data/pacbio/hifi/
INHIC=/nesi/nobackup/ga03048/data/HiC/weta/
DATAHIC=Weta_HiC_raw_all_  # suffix is _R1.fastq.gz, _R2.fastq.gz
OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/
OUTPRE=weta-hic-hifiasm
date=$(date)
##############

mkdir -p $OUTDIR
cd $OUTDIR

echo Starting hifiasm assembly for ${DATA} at $date

# Hi-C phasing with paired-end short reads in two FASTQ files
# -f38 applied to save mem on kmer counting
hifiasm --primary -t $SLURM_CPUS_PER_TASK -f38 -o ${OUTPRE}.asm \
--h1 ${INHIC}${DATAHIC}R1.fastq.gz --h2 ${INHIC}${DATAHIC}R2.fastq.gz \
${INDIR}weta_m64268e_210804_104249.hifi_reads.fastq.gz ${INDIR}weta_m64268e_210713_081012.hifi_reads.fastq.gz \
${INDIR}weta_m64267e_210804_084734.hifi_reads.fastq.gz ${INDIR}weta_m64267e_210802_234132.hifi_reads.fastq.gz \
${INDIR}weta_m64267e_210725_062308.hifi_reads.fastq.gz \
2> ${OUTPRE}.log

echo Finished at $date

awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.asm.hic.p_ctg.gfa > ${OUTPRE}.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.asm.hic.a_ctg.gfa > ${OUTPRE}.a_ctg.fa

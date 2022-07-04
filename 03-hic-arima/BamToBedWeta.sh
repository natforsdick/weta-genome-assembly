#!/bin/bash

################################
# Created 2020-11-26 by Nat Forsdick
# For preparing HiC data for input to SALSA
################################

module purge
module load BEDTools/2.29.2-GCC-9.2.0

LABEL='Weta_HiCmapped'
REP_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/weta/04_dedup/'
REF_DIR='/nesi/project/ga03048/results/redbean/weta-asm1-cns/'
REF='weta-asm1-wtdbg.cns.fa'

echo “Converting ${LABEL}\_rep1.bam to bed file”
bamToBed -i ${REP_DIR}${LABEL}\_rep1.bam > ${REP_DIR}${LABEL}\_rep1.bed
echo “converted”
echo “Sorting ${LABEL}\_rep1.bed”
sort -k 4 ${REP_DIR}${LABEL}\_rep1.bed > /nesi/nobackup/ga03048/tmp/tmp_bed${LABEL} && mv /nesi/nobackup/ga03048/tmp/tmp_bed${LABEL} ${REP_DIR}${LABEL}\_rep1.bed
echo “sorted”

# Need the .fai index file of the assembly - SALSA needs this to get the contig lengths.
if [ ! -e ${REF_DIR}${REF}.fai ]; then
module purge
module load SAMtools/1.10-GCC-9.2.0
echo “Making FAI from reference assembly”
cd $REF_DIR
samtools faidx ${REF_DIR}${REF}
echo “made”
else 
echo "Index file found"
fi


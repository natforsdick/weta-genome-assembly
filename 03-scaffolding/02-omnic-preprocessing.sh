#!/bin/bash -e

ml purge 
ml SAMtools/1.15.1-GCC-11.3.0 BWA/0.7.17-GCC-11.3.0

REF=02-weta-hic-hifiasm-p_ctg-purged # prefix
OUTDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/02-purge-dups/
TMPDIR=/nesi/nobackup/ga03048/tmp-omnic

mkdir $TMPDIR

cd $OUTDIR

samtools faidx $REF.fa

cut -f1,2 $REF.fa.fai > $REF.genome

bwa index $REF.fa

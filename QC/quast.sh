#!/bin/bash -e

################### 
# QUAST - assembly quality metrics
# Nat Forsdick, 2020-11-25
################### 

#### MODULES ####  
module purge
module load QUAST/5.0.2-gimkl-2018b
################### 

#### ENVIRONMENT #### 
datadir=/nesi/nobackup/ga03048/assemblies/SALSA/weta/05_SalsaScaff/
refdir=/nesi/project/ga03048/results/redbean/weta-asm1-cns/
ref=weta-asm1-wtdbg.cns.fa
outdir=/nesi/nobackup/ga03048/assemblies/SALSA/weta/05_SalsaScaff/
asmname=weta-asm1-wtdbg.cns.HiC_SALSA
################### 

#execute
if [ ! -e $outdir/$asmname ]; then
mkdir $outdir/$asmname
fi

cd $datadir
# Be sure to check that the correct .fa/.fasta suffix is used here.
quast.py -o $outdir/$asmname -t 8 --eukaryote ${refdir}${ref} ${datadir}${asmname}.fasta

#!/bin/bash -e

# Purge_dups pipeline
# Based on scripts by Sarah Bailey (UoA) to run the *purge_dups* pipeline (https://github.com/dfguan/purge_dups)
# Modified by Nat Forsdick, 2021-08-24

# step 01: align HiFi sequencing data to the assembly and generate a paf file
# Takes one parameter - PRI or ALT

#########
# MODULES
module purge
module load minimap2/2.20-GCC-9.2.0 
#########

#########
# PARAMS
INDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/01-assembly/
OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/02-purge-dups/
DATA=/nesi/project/ga03048/data/pacbio/hifi/
PRE=weta-hic-hifiasm # PREFIX - testing with single contig
PRI=p_ctg
ALT=a_ctg
R1=01- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
R2=02-

#########

if [ ! -e $OUTDIR ]; then
    mkdir -p $OUTDIR
    fi

cd $OUTDIR

if [ "$1" == "PRI" ]; then
# Testing first with just one fastq file
  minimap2 -x map-hifi -t $SLURM_CPUS_PER_TASK ${INDIR}${PRE}.${PRI}.fa ${DATA}*.fastq.gz | gzip -c - > ${R1}${PRE}-${PRI}-mapped.paf.gz
 
elif [ "$1" == "ALT" ]; then
  minimap2 -x map-hifi -t $SLURM_CPUS_PER_TASK ${OUTDIR}${R2}${PRE}.${ALT}.hap-merged.fa ${DATA}*.fastq.gz | gzip -c - > ${R1}${PRE}-${ALT}-merged-mapped.paf.gz

else
  minimap2 -x map-hifi -t $SLURM_CPUS_PER_TASK ${INDIR}${PRE}.fasta ${DATA}*.fastq.gz | gzip -c - > ${R1}${PRE}-mapped.paf.gz

fi

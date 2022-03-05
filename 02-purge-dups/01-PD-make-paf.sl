#!/bin/bash -e

#SBATCH --job-name=make-paf
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=22G
#SBATCH --ntasks=1
#SBATCH --profile=task 
#SBATCH --account=ga03048
#SBATCH --cpus-per-task=62

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
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
INDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/
OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/purge-dups/
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

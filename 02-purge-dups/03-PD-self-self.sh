#!/bin/bash -e

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# step 03: do a self-self alignment
# Takes one argument: PRI or ALT
#########
# MODULES
module purge
module load minimap2/2.20-GCC-9.2.0 
#########

#########
# PARAMS
OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/02-purge-dups/
PRE=weta-hic-hifiasm # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
R2=02-
#########

cd $OUTDIR
# -x asm5: intra-specific asm-to-asm alignment
if [ "$1" == "PRI" ]; then 
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}-${PRI}.split ${R1}${PRE}-${PRI}.split | gzip -c - > ${R1}${PRE}-${PRI}.split.self.paf.gz
elif [ "$1" == "ALT" ]; then
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}-${ALT}.split ${R1}${PRE}-${ALT}.split | gzip -c - > ${R1}${PRE}-${ALT}.split.self.paf.gz
else
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}.split ${R1}${PRE}.split | gzip -c - > ${R1}${PRE}.split.self.paf.gz
fi 

#!/bin/bash -e

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# step 06: modify cutoffs

##########
# PARAMS
PURGE_DUPS=/nesi/nobackup/ga03186/purge_dups/bin/
OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/purge-dups/
PRE=weta-hic-hifiasm- # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01-
R2=02- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
CUTOFFS="-l5 -m15 -u93"
##########

cd ${OUTDIR}
echo $CUTOFFS
${PURGE_DUPS}calcuts ${CUTOFFS} ${R1}${PRE}${PRI}-PB.stat > ${R2}${PRE}${PRI}-cutoffs

# Following this, you need to run steps 04-06 with $ROUND modified for new cutoffs.

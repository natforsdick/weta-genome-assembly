#!/bin/bash -e
# Generating FASTA from manually curated hic scaffolded assembly

##########
# PARAMS #
JUICER=/nesi/project/ga03186/scripts/Hi-C_scripts/yahs/juicer
REFDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/02-purge-dups/
REF=02-weta-hic-hifiasm-p_ctg-purged.fa
OUTDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/
CURATEDPREFIX=out_JBAT-2023-12-12
INPREFIX=out_JBAT

cd $OUTDIR
# -o = output prefix, followed by input review.assembly, input liftover.agp, reference
$JUICER post -o ${CURATEDPREFIX}-curated ${CURATEDPREFIX}.review.assembly ${INPREFIX}.liftover.agp ${REFDIR}${REF}

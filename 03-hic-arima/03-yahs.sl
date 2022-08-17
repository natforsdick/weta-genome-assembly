#!/bin/bash
#SBATCH --account=ga03048
#SBATCH --job-name=weta-yahs # job name (shows up in the queue)
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:20:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

################################
# Created 2020-11-26 by Nat Forsdick
# Passing aligned HiC Weta data to YAHS scaffolding genome assemblies
################################

# make output directory prior to running.
YAHS='/nesi/project/ga03186/scripts/Hi-C_scripts/yahs/yahs'
REF_DIR='/nesi/nobackup/ga03048/assemblies/hifiasm/01-assembly/'
REF='weta-hic-hifiasm.p_ctg.fa'
IN_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/hifiasm/04_dedup/'
IN_BAM='Weta-HiC-mapped_rep1.bed'
OUT_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/hifiasm/05-scaffolding/'

if [ ! -e ${OUT_DIR} ]; then
	mkdir -p ${OUT_DIR}
else
	echo "Found ${OUT_DIR}"
fi
cd ${OUT_DIR}

echo "Starting YAHS for ${IN_BAM} to scaffold ${REF}"
date

$YAHS ${REF_DIR}${REF} ${IN_DIR}${IN_BAM} -o yahs

echo "Completed YAHS scaffolding"
date

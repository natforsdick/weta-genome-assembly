#!/bin/bash
#SBATCH --account=ga03048
#SBATCH --job-name=weta-yahs # job name (shows up in the queue)
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

################################
# Created 2023-09-12 by Nat Forsdick
# Passing aligned HiC Weta data to YAHS scaffolding genome assemblies
################################

# make output directory prior to running.
YAHS='/nesi/project/ga03186/scripts/Hi-C_scripts/yahs/yahs'
REF_DIR='/nesi/nobackup/ga03048/weta/assemblies/hifiasm/02-purge-dups/'
REF='02-weta-hic-hifiasm-p_ctg-purged.fa'
IN_DIR='/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/01-mapping/'
IN_BAM='02-weta-hic-hifiasm-p_ctg-purged-AG1149-mapped.PT.bam'
OUT_DIR='/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/'

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

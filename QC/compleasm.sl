#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J compleasm # job name (shows up in the queue)
#SBATCH -c 28
#SBATCH --mem=58GB #
#SBATCH --time=02:00:00 #Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out #
#SBATCH --error %x.%j.err #

# compleasm - 2023-06-15
# Genome assembly QC - BUSCO alternative

# PARAMS
INDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/omnic-r2/yahs/
OUTDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/BUSCO/
DB=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/BUSCO/insecta_odb12
ASM=out_JBAT-2023-12-12-curated-AG1149-mapped.PT-yahs_scaffolds_final.fa

asm=$(basename $ASM .fa)

cd $INDIR
module purge
module load compleasm/0.2.7-gimkl-2022a

compleasm run -a ${INDIR}${ASM} -o ${OUTDIR}compleasm-${asm} -t 24 -l insecta -L $DB

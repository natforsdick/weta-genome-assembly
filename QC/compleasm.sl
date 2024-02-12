#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J compleasm # job name (shows up in the queue)
#SBATCH -c 32
#SBATCH --mem=64GB #
#SBATCH --time=01:30:00 #Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out #
#SBATCH --error %x.%j.err #

# compleasm - 2023-06-15
# Genome assembly QC - BUSCO alternative

# PARAMS
INDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/
OUTDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/BUSCO/
DB=/nesi/project/ga03048/insecta_odb10/insecta_odb10
ASM=out_JBAT-2023-12-12-curated.FINAL.fa

asm=$(basename $ASM .fa)

cd $INDIR
module purge
module load compleasm/0.2.2-gimkl-2022a

compleasm.py run -a ${INDIR}${ASM} -o ${OUTDIR}compleasm-${asm} -t 18 -l insecta -L $DB

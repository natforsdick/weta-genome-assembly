#!/bin/bash -e

#SBATCH --job-name=minimap-aln
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=50G #Used 26 for whole genomes
#SBATCH --profile=task 
#SBATCH --account=ga03048
#SBATCH --cpus-per-task=40 # Used 32 for whole genomes

# align-genomes.sl
# Align genomes to one another using minimap
# Nat Forsdick, 2021-09-01

#########
# MODULES
module purge
module load minimap2/2.20-GCC-9.2.0 
#########

INDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/02-purge-dups/
QUERY=01-weta-hic-hifiasm.a_ctg.hap-merged
REFDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/02-purge-dups/
REF=01-weta-hic-hifiasm-p_ctg-purged
OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/04-alignment/

cd $OUTDIR

echo "Aligning $QUERY against $REF reference genome"

# To index reference genome the first time you run this - can then just call the index ref.mmi following this
if [ ! -e ${REFDIR}${REF}.mmi ]; then
	echo "Index file of reference does not exist: creating index"
	minimap2 -t $SLURM_CPUS_PER_TASK -d ${REFDIR}${REF}.mmi ${REFDIR}${REF}.fa
fi
echo "Aligning $QUERY to $REF"
minimap2 -ax asm5 -t 32 ${REFDIR}${REF}.mmi ${INDIR}${QUERY}.fa > pri-pur-v-alt-hap-merged.paf

echo "Collecting stats"
$HOME/bin/k8 /nesi/project/ga03186/HiFi-scripts/paftools.js stat pri-pur-v-alt-hap-merged.paf > pri-pur-v-alt-hap-merged-stat.out

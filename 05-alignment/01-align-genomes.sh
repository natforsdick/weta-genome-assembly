#!/bin/bash -e

# align-genomes.sl
# Align genomes to one another using minimap
# Nat Forsdick, 2021-09-01

#########
# MODULES
module purge
module load minimap2/2.24-GCC-11.3.0
#########

INDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/omnic-r2/yahs/
QUERY=out_JBAT-2023-12-12-curated-AG1149-mapped.PT-yahs_scaffolds_final
REFDIR=/nesi/nobackup/ga03048/weta/reference-genomes/
REF=GCA_946902985.2_iqMecThal1.2_genomic
OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/04-alignment/

export TMPDIR=/nesi/nobackup/ga03048/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR=$TMPDIR

cd $OUTDIR

echo "Aligning $QUERY against $REF reference genome"

# To index reference genome the first time you run this - can then just call the index ref.mmi following this
if [ ! -e ${REFDIR}${REF}.mmi ]; then
	echo "Index file of reference does not exist: creating index"
	minimap2 -t $SLURM_CPUS_PER_TASK -d ${REFDIR}${REF}.mmi ${REFDIR}${REF}.fna
fi
echo "Aligning $QUERY to $REF"
minimap2 -ax asm5 -t 32 ${REFDIR}${REF}.mmi ${INDIR}${QUERY}.fa > GCA_946902985.2-curated-AG1149-mapped.PT-yahs_scaffolds.paf

echo "Collecting stats"
paftools.js stat GCA_946902985.2-curated-AG1149-mapped.PT-yahs_scaffolds.paf > GCA_946902985.2-curated-AG1149-mapped.PT-yahs_scaffolds-stats.out

ml purge && ml miniasm/0.3-20191007-GCC-11.3.0

# -m min match length, -w image width
minidot -m 5000 -w 1000 GCA_946902985.2-curated-AG1149-mapped.PT-yahs_scaffolds.paf

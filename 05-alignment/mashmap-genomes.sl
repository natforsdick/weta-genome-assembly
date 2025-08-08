#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J mashmap
#SBATCH --time 01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err

# mashmap-genomes.sl
# Nat Forsdick, 2024-02-14
# aligning assemblies to reference genomes

REFDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/omnic-r2/yahs/
REF=out_JBAT-2023-12-12-curated-AG1149-mapped.PT-yahs_scaffolds_final.fa
QUERYDIR=/nesi/nobackup/ga03048/weta/reference-genomes/
QUERY=GCA_946902985.2_iqMecThal1.2_genomic.fna
OUTFILE=curated-AG1149-mapped.PT-yahs_scaffolds_final-v-GCA_946902985.2
cd $QUERYDIR

ml purge
ml MashMap/3.0.4-Miniconda3

mashmap -r ${REFDIR}${REF} -q ${QUERYDIR}${QUERY} --perc_identity 80 -f one-to-one -t 10 -o ${OUTFILE}-mashmap
generateDotPlot png medium ${OUTFILE}-mashmap.out

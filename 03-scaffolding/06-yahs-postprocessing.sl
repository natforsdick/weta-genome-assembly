#!/bin/bash
#SBATCH --account=ga03048
#SBATCH --job-name=make-juice # job name (shows up in the queue)
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=4:00:00 #Walltime (HH:MM:SS) # 8 hrs to run whole two step pipeline
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err 

# Generating .hic file for visualisation and manual curation

##########
# PARAMS #
YAHSJUICE=/nesi/project/ga03186/scripts/Hi-C_scripts/yahs/juicer
OUTDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/omnic-r2/yahs/
JUICER=/nesi/nobackup/ga03048/juicer/scripts/juicer_tools.1.9.9_jcuda.0.8.jar
REF_DIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/
REF=out_JBAT-2023-12-12-curated.FINAL.fa
REFPRE=02-weta-hic-hifiasm-p_ctg-purged
TMPDIR=/nesi/nobackup/ga03048/tmp
##########

cd $OUTDIR
export TMPDIR

#echo 'generating hic contact map'
#$YAHSJUICE pre -a -o out_JBAT_omnic-r2 yahs.bin out_JBAT-2023-12-12-curated-AG1149-mapped.PT-yahs_scaffolds_final.agp \
#${REF_DIR}${REF}.fai > out_JBAT_omnic-r2.log 2>&1
#echo done step 1

echo 'running juicer_tools pre'
(java -jar -Xmx48G $JUICER pre out_JBAT_omnic-r2.txt out_JBAT_omnic-r2.hic.part <(cat out_JBAT_omnic-r2.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT_omnic-r2.hic.part out_JBAT_omnic-r2.hic)

#mv out_JBAT.hic.part out_JBAT.hic

echo done step 2



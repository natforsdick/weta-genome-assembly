#!/bin/bash -e
#SBATCH --account=ga03048
#SBATCH --job-name=omnic-map # job name (shows up in the queue)
#SBATCH --cpus-per-task=64 # mapping can use 18, subsequent processing requires 6
#SBATCH --mem=35G
#SBATCH --time=2-10:00:00 #3-00:00:00 #Walltime (HH:MM:SS) # Total processing minimum 12 hrs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

ml purge
ml SAMtools/1.15.1-GCC-11.3.0 BWA/0.7.17-GCC-11.3.0 

#########
# PARAMS
PREFIX=out_JBAT-2023-12-12-curated-AG1149
INDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/2023-OmniC-out/
OUTDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/omnic-r2/
REFDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/02-scaffolding/
REFPRE=out_JBAT-2023-12-12-curated.FINAL
REF=${REFDIR}${REFPRE}.fa
TMPDIR="/nesi/nobackup/ga03186/tmp-omnic-${SLURM_JOB_ID}"
FQ1=_R1_001_clean.fastq.gz
FQ2=_R2_001_clean.fastq.gz
CPU=32
########

mkdir $TMPDIR
cd $INDIR

# alignment -T0 = all alignments to generate stats, quality filtering later
# we have to process all files together
echo aligning $file
bwa mem -5SP -T0 -t $CPU $REF <(zcat AG1149-01_S1_L001_R1_001_clean.fastq.gz AG1149-01_S1_L002_R1_001_clean.fastq.gz AG1149-02_S2_L001_R1_001_clean.fastq.gz AG1149-02_S2_L002_R1_001_clean.fastq.gz AG1149-03_S3_L001_R1_001_clean.fastq.gz AG1149-03_S3_L002_R1_001_clean.fastq.gz) \
<(zcat AG1149-01_S1_L001_R2_001_clean.fastq.gz AG1149-01_S1_L002_R2_001_clean.fastq.gz AG1149-02_S2_L001_R2_001_clean.fastq.gz AG1149-02_S2_L002_R2_001_clean.fastq.gz AG1149-03_S3_L001_R2_001_clean.fastq.gz AG1149-03_S3_L002_R2_001_clean.fastq.gz) \
-o ${OUTDIR}${PREFIX}-aligned.sam
echo aligment complete

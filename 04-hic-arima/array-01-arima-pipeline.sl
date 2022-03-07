#!/bin/bash 
#SBATCH --account=ga03048
#SBATCH --job-name=01-arima-indexing # job name (shows up in the queue)
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=01:30:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.%a.out # CHANGE number for new run
#SBATCH --error %x.%j.%a.err #  CHANGE number for new run
#SBATCH --array=1-2
#SBATCH --profile=task
################################


##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
# https://github.com/ArimaGenomics/mapping_pipeline/ #
##############################################

# This is the ARIMA genomics pipeline, modified by Nat Forsdick for use with wētā data.

# Below find the commands used to map HiC data.

# Replace the variables at the top with the correct paths for the locations of files/programs on your system.

# This bash script will map one paired end HiC dataset (read1 & read2 fastqs). 

# Modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

############
# MODULES
############
module purge
module load BWA/0.7.17-GCC-9.2.0 picard/2.21.8-Java-11.0.4 SAMtools/1.13-GCC-9.2.0 samblaster/0.1.26-GCC-9.2.0
############

############
# PARAMS
############
HIC='Weta_HiC_raw_all_R' #'basename_of_fastq_files'
LABEL='Weta-HiC-mapped' #'overall_exp_name'
BWA='bwa' #'/software/bwa/bwa-0.7.12/bwa'
SAMTOOLS='samtools' #'/software/samtools/samtools-1.3.1/samtools'
IN_DIR='/nesi/nobackup/ga03048/data/HiC/weta/' # '/path/to/gzipped/fastq/files'
REF_DIR='/nesi/nobackup/ga03048/assemblies/hifiasm/01-assembly/'
REF='/nesi/nobackup/ga03048/assemblies/hifiasm/01-assembly/weta-hic-hifiasm.p_ctg.fa' #'/path/to/reference_sequences/reference_sequences.fa'
FAIDX='$REF.fai'
PREFIX='weta-hic-hifiasm.p_ctg' #'bwa_index_name'
RAW_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/hifiasm-raw/01_mapped/' #'/path/to/write/out/bams'
FILT_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/hifiasm-raw/02_filtered/' #'/path/to/write/out/filtered/bams'
FILTER='/nesi/project/ga03048/scripts/genome-assembly/04-hic-arima/filter_five_end.pl' #'/path/to/filter_five_end.pl'
COMBINER='/nesi/project/ga03048/scripts/genome-assembly/04-hic-arima/two_read_bam_combiner.pl' #'/path/to/two_read_bam_combiner.pl'
STATS='/nesi/project/ga03048/scripts/genome-assembly/04-hic-arima/get_stats.pl' #'/path/to/get_stats.pl'
PICARD='/opt/nesi/mahuika/picard/2.21.8-Java-11.0.4/picard.jar'
TMP_DIR='/nesi/nobackup/ga03048/tmp/' #'/path/to/write/out/temporary/files'
PAIR_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/hifiasm-raw/03_paired/' #'/path/to/write/out/paired/bams'
REP_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/hifiasm-raw/04_dedup/' #'/path/to/where/you/want/deduplicated/files'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/hifiasm-raw/05_merged/' #'/path/to/final/merged/alignments/from/any/biological/replicates'
MAPQ_FILTER=10
CPU=$SLURM_CPUS_PER_TASK
############

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
cd $REF_DIR
if [ -f ${REF}.amb ]; then
    echo "${REF} index found"
    else
	bwa index -a bwtsw -p $PREFIX $REF
	echo "Finished indexing $REF"
fi

echo "### Step 1.A: FASTQ to BAM (1st)"

samplesheet="${IN_DIR}samplesheet.txt"

echo "r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`"

r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`

echo "bwa mem -t $CPU -5SP $REF ${IN_DIR}${r1} ${IN_DIR}${r2} | samblaster |\
 samtools view -@ $CPU -buSh -F 2316 - > ${RAW_DIR}${r1}.bam"

bwa mem -t $CPU -5SP $REF ${IN_DIR}${r1} ${IN_DIR}${r2} | samblaster |\
 samtools view -@ $CPU -buSh -F 2316 - > ${RAW_DIR}${r1}.bam

echo "Finished step 1.A"



#!/bin/bash 
#SBATCH --account=ga03048
#SBATCH --job-name=arima_mapping # job name (shows up in the queue)
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=20:00:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

################################


######################################################
#      ARIMA GENOMICS MAPPING PIPELINE 02/08/2019    #
# https://github.com/ArimaGenomics/mapping_pipeline/ #
######################################################
# Modified by Natalie Forsdick 2020-11-24
# For use with Weta HiC data
######################################################
#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

# Modules
module purge
module load BWA/0.7.17-gimkl-2017a picard/2.21.8-Java-11.0.4 SAMtools/1.10-GCC-9.2.0 


SRA='Weta_HiC_raw_all_R' #'basename_of_fastq_files'
LABEL='Weta_HiCmapped' #'overall_exp_name'
BWA='bwa' #'/software/bwa/bwa-0.7.12/bwa'
SAMTOOLS='samtools' #'/software/samtools/samtools-1.3.1/samtools'
IN_DIR='/nesi/nobackup/ga03048/data/HiC/weta' # '/path/to/gzipped/fastq/files'
REF='/nesi/project/ga03048/results/redbean/weta-asm1-cns/weta-asm1-wtdbg.cns.fa' #'/path/to/reference_sequences/reference_sequeneces.fa'
FAIDX='$REF.fai'
PREFIX='weta-asm1-wtdbg.cns' #'bwa_index_name'
RAW_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/weta/01_mapped' #'/path/to/write/out/bams'
FILT_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/weta/02_filtered' #'/path/to/write/out/filtered/bams'
FILTER='/nesi/project/ga03048/scripts/arima_pipeline/filter_five_end.pl' #'/path/to/filter_five_end.pl'
COMBINER='/nesi/project/ga03048/scripts/arima_pipeline/two_read_bam_combiner.pl' #'/path/to/two_read_bam_combiner.pl'
STATS='/nesi/project/ga03048/scripts/arima_pipeline/get_stats.pl' #'/path/to/get_stats.pl'
PICARD='/opt/nesi/mahuika/picard/2.21.8-Java-11.0.4/picard.jar'
TMP_DIR='/nesi/nobackup/ga03048/tmp' #'/path/to/write/out/temporary/files'
PAIR_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/weta/03_paired' #'/path/to/write/out/paired/bams'
REP_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/weta/04_dedup' #'/path/to/where/you/want/deduplicated/files'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='/nesi/nobackup/ga03048/assemblies/SALSA/weta/05_merged' #'/path/to/final/merged/alignments/from/any/biological/replicates'
MAPQ_FILTER=10
CPU=16

#echo "### Step 0: Check output directories exist & create them as needed"
#[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
#[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
#[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
#[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
#[ -d $REP_DIR ] || mkdir -p $REP_DIR
#[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

#echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
#$BWA index -a bwtsw -p $PREFIX $REF

#echo "### Step 1.A: FASTQ to BAM (1st)"
#$BWA mem -t $CPU $REF ${IN_DIR}/${SRA}1.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\1.bam
#echo " Ready for step 1.A: $BWA mem -t $CPU $REF ${IN_DIR}/${SRA}\1.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > ${RAW_DIR}/${SRA}\1.bam " 
#echo "### Step 1.B: FASTQ to BAM (2nd)"
#$BWA mem -t $CPU $REF ${IN_DIR}/${SRA}\1.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > ${RAW_DIR}/${SRA}\2.bam

#echo "### Step 2.A: Filter 5' end (1st)"
#samtools view -@ $CPU -h ${RAW_DIR}/${SRA}1.bam | perl $FILTER | samtools view -@ $CPU -Sb - > ${FILT_DIR}/${SRA}1.bam

#echo "### Step 2.B: Filter 5' end (2nd)"

#samtools view -@ $CPU -h ${RAW_DIR}/${SRA}2.bam | perl $FILTER | samtools view -@ $CPU -Sb - > ${FILT_DIR}/${SRA}2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER ${FILT_DIR}/${SRA}1.bam ${FILT_DIR}/${SRA}2.bam samtools $MAPQ_FILTER | samtools view -@ $CPU -bS -t $FAIDX - | samtools sort -@ $CPU -o ${TMP_DIR}/${SRA}.bam -

echo "### Step 3.B: Add read group"
java -Xmx20G -Djava.io.tmpdir=${TMP_DIR} -jar $PICARD AddOrReplaceReadGroups INPUT=${TMP_DIR}/${SRA}.bam OUTPUT=${PAIR_DIR}/${SRA}.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none
#echo "Completed Step 3.B"
###############################################################################################################################################################
###                                           How to Accommodate Technical Replicates                                                                       ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                                                    ###
### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
### The code below is an example of how to merge technical replicates.                                                                                      ###
###############################################################################################################################################################
#	REP_NUM=X #number of the technical replicate set e.g. 1
#	REP_LABEL=$LABEL\_rep$REP_NUM
#	INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as technical replicates
#   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
java -Xmx20G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=${TMP_DIR} -jar $PICARD MarkDuplicates INPUT=${PAIR_DIR}/${SRA}.bam OUTPUT=${REP_DIR}/${REP_LABEL}.bam METRICS_FILE=${REP_DIR}/metrics.${REP_LABEL}.txt TMP_DIR=${TMP_DIR} ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index -@ $CPU ${REP_DIR}/${REP_LABEL}.bam

perl $STATS ${REP_DIR}/${REP_LABEL}.bam > ${REP_DIR}/${REP_LABEL}.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

#########################################################################################################################################
###                                       How to Accommodate Biological Replicates                                                    ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                              ###
### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
### The code below is an example of how to merge biological replicates.                                                               ###
#########################################################################################################################################
#
#	INPUTS_BIOLOGICAL_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as biological replicates
#   example bash array - INPUTS_BIOLOGICAL_REPS=('INPUT=A_rep1.bam' 'INPUT=A_rep2.bam' 'INPUT=A_rep3.bam')
#
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$MERGE_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
#
#	$SAMTOOLS index $MERGE_DIR/$LABEL.bam

# perl $STATS $MERGE_DIR/$LABEL.bam > $MERGE_DIR/$LABEL.bam.stats

# echo "Finished Mapping Pipeline through merging Biological Replicates"

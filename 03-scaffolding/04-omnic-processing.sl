#!/bin/bash -e
#SBATCH --account=ga03048
#SBATCH --job-name=omnic-map # job name (shows up in the queue)
#SBATCH --cpus-per-task=12 # mapping can use 18, subsequent processing requires 6
#SBATCH --mem=10G ## mod based on input size
#SBATCH --time=6-00:00:00 #3-00:00:00 #Walltime (HH:MM:SS) # Total processing minimum 12 hrs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run
#SBATCH --profile=task=300

ml purge 
ml SAMtools/1.15.1-GCC-11.3.0 BWA/0.7.17-GCC-11.3.0 pairtools/1.0.2-gimkl-2022a-Python-3.10.5

#########
# PARAMS
PREFIX=02-weta-hic-hifiasm-p_ctg-purged-AG1149
INDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/03-scaffolding/2023-OmniC/01-mapping/
REFDIR=/nesi/nobackup/ga03048/weta/assemblies/hifiasm/02-purge-dups/
REFPRE=02-weta-hic-hifiasm-p_ctg-purged
REF=${REFDIR}${REFPRE}.fa
TMPDIR="/nesi/nobackup/ga03048/tmp-omnic-${SLURM_JOB_ID}"
########

mkdir $TMPDIR
cd $INDIR

# find ligation junctions
echo finding ligation junctions
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 24 \
--nproc-out 24 --chroms-path ${REFDIR}${REFPRE}.genome ${PREFIX}-aligned.sam > ${PREFIX}-parsed.pairsam
echo ligation junctions found

# sort pairsam
echo sorting pairsam
pairtools sort --nproc 24 --tmpdir=$TMPDIR ${PREFIX}-parsed.pairsam > ${PREFIX}-sorted.pairsam

# remove duplicates
echo removing duplicates
pairtools dedup --nproc-in 24 --nproc-out 24 --mark-dups --output-stats ${PREFIX}-stats.txt \
--output ${PREFIX}-dedup.pairsam ${PREFIX}-sorted.pairsam

# split .bam, .pairs
echo splitting bam
pairtools split --nproc-in 24 --nproc-out 24 --output-pairs ${PREFIX}-mapped.pairs \
--output-sam ${PREFIX}-unsorted.bam ${PREFIX}-dedup.pairsam

# sort bam
echo sorting bam
samtools sort -@24 -T ${TMPDIR}tempfile.bam -o ${PREFIX}-mapped.PT.bam ${PREFIX}-unsorted.bam

# index bam
echo indexing final bam
samtools index ${PREFIX}-mapped.PT.bam

if [ -f ${PREFIX}-mapped.PT.bam ]
then
echo "pipeline completed"
else
echo "pipeline not complete"
fi

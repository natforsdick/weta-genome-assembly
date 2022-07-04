#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J weta-polish
#SBATCH --cpus-per-task=72 # full run requires 72
#SBATCH --mem=20G # full run requires 75G
#SBATCH --time=04:00:00 # full run requires 10 hrs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

#######
# PARAMS
asmdir=/nesi/nobackup/ga03048/assemblies/hifiasm/purge-dups/
fo=02-weta-hic-hifiasm-p_ctg-purged # testing with subset longest-contig.weta-hic-hifiasm
datadir=/nesi/project/ga03048/data/pacbio/hifi/
infastq=weta_m64267e_210804_084734.hifi_reads #test
outdir=/nesi/nobackup/ga03048/assemblies/hifiasm/polishing/
######

ml purge
ml wtdbg/2.5-GCC-9.2.0 BWA/0.7.17-GCC-9.2.0 \
SAMtools/1.13-GCC-9.2.0 minimap2/2.20-GCC-9.2.0

if [ ! -e ${outdir} ]; then
	mkdir -p $outdir
fi

cd $outdir

date
echo "minimap2 -I 24G -t 60 -ax map-hifi -r2k ${asmdir}${fo}.fa ${datadir}*.fastq.gz |\
samtools sort -@4 -o ${fo}.bam"
#minimap2 -I 24G -t 60 -ax map-hifi -r2k ${asmdir}${fo}.fa ${datadir}*.fastq.gz |\
#samtools sort -@4 -o ${fo}.bam

date
echo "samtools view -F0x900 ${fo}.bam | wtpoa-cns -t 60 \
-d ${asmdir}${fo}.fa -i - -fo ${fo}.cns.fa"
samtools view -F0x900 ${fo}.bam | wtpoa-cns -t 60 \
-d ${asmdir}${fo}.fa -i - -fo ${fo}.cns.fa

date
echo "bwa index ${fo}.cns.fa"
bwa index ${fo}.cns.fa

# If you wanted to do this with short-read data (probably not recommended with HiFi, as you
# may incur mismapping across repetitive regions from Illumina)≈
#bwa mem -t $SLURM_CPUS_PER_TASK ${fo}.cns.fa ${srdir}/*fastq.gz | \
#samtools sort -O SAM | wtpoa-cns -t 24-x sam-sr \
#-d ${fo}.cns.fa -i - -fo ${fo}.srp.fa 
# wtpoa-cns calls in short-read data for polishing - would want to polish with CCS first



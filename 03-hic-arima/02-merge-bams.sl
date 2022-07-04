#!/bin/bash 
#SBATCH --account=ga03048
#SBATCH --job-name=02-merge # job name (shows up in the queue)
#SBATCH --cpus-per-task=16
#SBATCH --mem=2G
#SBATCH --time=01:00:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run
#SBATCH --profile=task
################################

ml SAMtools/1.13-GCC-9.2.0

cd /nesi/nobackup/ga03048/assemblies/SALSA/hifiasm-raw/01_mapped/

samtools merge -@ $SLURM_CPUS_PER_TASK Weta_HiC_raw_all.bam Weta_HiC_raw_all_R1.fq.*.bam
echo done

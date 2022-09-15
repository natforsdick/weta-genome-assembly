#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J RagTag
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=22G
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err
#SBATCH --profile=task

# RagTagPatch.sl
# Nat Forsdick, 2022-09-15
# Joining and filling a draft assembly using a reference assembly
# In this case, weta redbean to fill the HiFiasm output

# PARAMS #
QUERY=/nesi/project/ga03048/results/redbean/weta-asm1-cns/weta-asm1-wtdbg.cns.fa
TARGET=/nesi/nobackup/ga03048/assemblies/hifiasm/02-purge-dups/01-weta-hic-hifiasm-p_ctg-purged.fa
OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/06-GapClose/

# ENVIRONMENT #
ml Miniconda3/4.12.0 minimap2/2.24-GCC-9.2.0
source /opt/nesi/CS400_centos7_bdw/Miniconda3/4.8.3/etc/profile.d/conda.sh
conda activate ragtag

cd $OUTDIR
ragtag.py patch $TARGET $QUERY -o $OUTDIR -t 12 --aligner minimap2
conda deactivate
ml purge

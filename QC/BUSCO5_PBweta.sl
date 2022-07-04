#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J BUSCO5_altpur # job name (shows up in the queue)
#SBATCH -c 32
#SBATCH --mem=38GB
#SBATCH --partition=large
#SBATCH --time=14:00:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # 
#SBATCH --error %x.%j.err #  


###################
# BUSCO - INSECTS - 2020-12-13
# Nat Forsdick
###################

# Load modules
module purge
module load BUSCO/5.2.2-gimkl-2020a
#cp -r $AUGUSTUS_CONFIG_PATH ./MyAugustusConfig

# Set up environment
#export AUGUSTUS_CONFIG_PATH=/nesi/project/ga03048/scripts/QC/MyAugustusConfig
#export BUSCO_CONFIG_FILE="/nesi/project/ga03048/scripts/config.ini"

OUTDIR=/nesi/nobackup/ga03048/assemblies/hifiasm/BUSCO/
IN_DIR=/nesi/nobackup/ga03048/assemblies/hifiasm/02-purge-dups/
samplist='01-weta-hic-hifiasm-a_ctg-purged' #'weta-hic-hifiasm.cns'
INSECT_DB=/nesi/project/ga03048/insecta_odb10

for samp in $samplist
do

mkdir -p $OUTDIR
cd $OUTDIR

echo "Starting BUSCO for ${samp}"
# -f = force, -r = restart
	busco -i ${IN_DIR}${samp}.fa -o BUSCO5_${samp} -f --offline -l ${INSECT_DB} -m geno -c $SLURM_CPUS_PER_TASK
	echo "Finished BUSCO for ${samp}"
done

# To make BUSCO plots:
# ml BUSCO/5.2.2-gimkl-2020a R/4.1.0-gimkl-2020a
# generate_plot.py -wd ./



#!/bin/bash -e
#SBATCH -A ga03048
#SBATCH -J hic_qc
#SBATCH -c 2
#SBATCH --mem=500M
#SBATCH --partition=large
#SBATCH --time=00:03:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out #
#SBATCH --error %x.%j.err #

###################
# Run HiC_QC
# Nat Forsdick, 2020-11-27
###################

# This script takes two arguments, $1 : the path to the input directory, /
# and $2 : the prefix of the input file.
# e.g. to execute:

bash ./hic-qc.sh /nesi/nobackup/ga03048/results/hic-qc-kaki/ Kaki_HiC_mapped

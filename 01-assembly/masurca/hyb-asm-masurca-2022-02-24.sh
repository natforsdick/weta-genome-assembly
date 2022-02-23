#!/bin/bash

#SBATCH --account ga03048
#SBATCH --job-name masurca-weta
#SBATCH --time 00:05:00
#SBATCH --mem 4G
#SBATCH --cpus-per-task 8
#SBATCH --qos=debug
#SBATCH --output hifiasm.%j.out
#SBATCH --error hifiasm.%j.err
#SBATCH --profile=task

# Test script to run masurca pipeline for weta

#load modules
module load MaSuRCA/4.0.5-gimkl-2020a

#create analysis env
analysisdir="/nesi/nobackup/ga03048/assemblies/masurca"
datadir=${analysisdir}/"00_rawdata"

# Move to analysis dir
cd ${analysisdir}

#execute
masurca nf_config.txt

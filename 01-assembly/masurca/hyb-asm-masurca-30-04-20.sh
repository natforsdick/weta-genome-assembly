#!/bin/bash

#SBATCH --account landcare02543
#SBATCH --job-name masurca-assemble-weta
#SBATCH --time 00:05:00
#SBATCH --mem 4G
#SBATCH --cpus-per-task 8
#SBATCH --qos=debug


# Test script to run canu  pipeline for weta

#load modules
module load MaSuRCA/3.3.4-gimkl-2018b

#create analysis env
analysisdir="/nesi/nobackup/landcare00070/checkpoint/manpreet.dhami/weta/01_masurca"
datadir=${analysisdir}/"00_rawdata"

# Move to analysis dir
cd ${analysisdir}

#execute
masurca sr_config.txt

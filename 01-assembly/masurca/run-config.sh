#!/bin/bash

# Make MaSuRCA config for weta

# MODULES
module load MaSuRCA/4.0.5-gimkl-2020a

# PARAMS
analysisdir="/nesi/nobackup/ga03048/assemblies/masurca"
datadir=${analysisdir}/"00_rawdata"

cd ${analysisdir}

masurca /nesi/project/ga03048/scripts/genome-assembly/01-assembly/masurca/nf_config.txt

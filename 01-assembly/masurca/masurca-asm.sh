#!/bin/bash -e

#############
# run-masurca.sl
# Nat Forsdick, 2021-11-01
# Masurca hybrid genome assembly
#############

############
# MODULES
ml purge
ml MaSuRCA/4.0.9-gimkl-2020a
ml minimap2/2.20-GCC-9.2.0
############

echo Beginnning MaSuRCA at
date

cd /nesi/nobackup/ga03048/assemblies/masurca4.0.9/

bash ./assemble.sh > 2022-04-26-assemble.out

echo Finishing MaSuRCA at
date

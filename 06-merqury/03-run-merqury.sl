#!/bin/bash -e 

#SBATCH -A ga03048
#SBATCH -J merqury
#SBATCH --mem=16G
#SBATCH -c 24
#SBATCH --time=00:40:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --profile=task

# run-merqury.sl
# N Forsdick, 2022-04-04
# Takes 3 params: 1) path to assembly directory, 2) primary assembly, 3) alternate assembly
##########
# PARAMS #
##########
genome=weta
outdir=/nesi/nobackup/ga03048/assemblies/hifiasm/05-merqury/
asmdir=$1
asm1=$2
asm2=$3
##########

cd $outdir

asmname=$(basename "$asm1")
asmname=${asmname%.*}

mkdir -p $asmname
cd $asmname

if [ -e "$asm1" ] && [ -e "$asm2" ]; then
	echo "${asm1} and ${asm2} exist"
else
ln -s ${asmdir}${asm1} ${asm1}
ln -s ${asmdir}${asm2} ${asm2}
fi

module purge
module load R/4.1.0-gimkl-2020a BEDTools/2.29.2-GCC-9.2.0 SAMtools/1.13-GCC-9.2.0

$MERQURY/merqury.sh ../${genome}.meryl ${asm1} ${asm2} ${asmname}-merqury

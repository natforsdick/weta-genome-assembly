#!/bin/bash -e

#SBATCH --job-name=self-self
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=46G
#SBATCH --ntasks=1
#SBATCH --profile=task
#SBATCH --account=ga03048
#SBATCH --cpus-per-task=14

# Takes one param: PRI/ALT

bash ./03-PD-self-self.sh PRI

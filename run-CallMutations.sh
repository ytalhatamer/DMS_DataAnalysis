#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=256GB
#SBATCH --error=jobs/job%j.err
#SBATCH --output=jobs/job%j.out
inputfolder=$1
## Count number of cores to distribute the jobs
numcores="$(grep -c ^processor /proc/cpuinfo)"
## One core is left idle for operating system
echo "Number of Cores used: " $numcores
echo "Input folder for analysis is: " $inputfolder

python PoolJob.py $numcores CallMutations $inputfolder

python PopulateDataCalculateEnrichment.py $inputfolder

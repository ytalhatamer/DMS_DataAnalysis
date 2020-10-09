#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=super
#SBATCH --error=jobs/job%j.err
#SBATCH --output=jobs/job%j.out
inputfolder=$1
## Count number of cores to distribute the jobs
numcores="$(grep -c ^processor /proc/cpuinfo)"
## One core is left idle for operating system
echo $numcores
echo $inputfolder
python PoolJob.py $numcores Flash $inputfolder

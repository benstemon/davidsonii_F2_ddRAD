#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_stacks


cd $SLURM_SUBMIT_DIR

jobscriptpath="/work/bs66/davidsonii_mapping/preprocessing_v3/stacks_jobscripts"

for i in $jobscriptpath/*.sh;
do
    sbatch $i
done

#!/bin/sh

#SBATCH -N 1
#SBATCH -n 10 
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_stacks


cd $SLURM_SUBMIT_DIR


module load stacks/gcc/2.41


infilepath="/work/bs66/davidsonii_mapping/preprocessing_v2/fastp_outfiles"
outfilepath="/work/bs66/davidsonii_mapping/preprocessing_v2/stacks_output"


for r1in in $infilepath/*R1.fq.gz;
do
    r2in="${r1in/_R1./_R2.}"
    process_radtags --paired -1 $r1in -2 $r2in -i gzfastq -o $outfilepath -c -q -w 0.15 -s 20 --len_limit 30 --disable_rad_check
done

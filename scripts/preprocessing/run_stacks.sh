#!/bin/sh

#SBATCH -N 1
#SBATCH -n 10 
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_stacks


cd $SLURM_SUBMIT_DIR


module load stacks/gcc/2.41


infilepath="/work/bs66/davidsonii_mapping/preprocessing_v1/fastp_outfiles"
outfilepath="/work/bs66/davidsonii_mapping/preprocessing_v1/stacks_output"


for r1in in $infilepath/PopF2_R1_*;
do
    r2in="${r1in/_R1_/_R2_}"
    process_radtags --paired -1 $r1in -2 $r2in -i gzfastq -o $outfilepath --renz_1 'ecoRI' --renz_2 'mspI' -c -q -r -w 0.15 -s 10 --len_limit 30
done

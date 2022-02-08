#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_fastp

cd $SLURM_SUBMIT_DIR


#path to fastp
fastpdir="/work/bs66/software"

#path to demultiplexed files
demuxdir="/work/bs66/davidsonii_mapping/demuxed"


#fastp for loop
for r1in in $demuxdir/*.F.fq.gz; 
do
    r2in="${r1in/F.fq.gz/R.fq.gz}"
    r1out="${r1in##*/}"
    r2out="${r1out/F.fq.gz/R.fq.gz}"
    $fastpdir/./fastp -i "$r1in" -I "$r2in" -o "${r1out/.F./_R1.}" -O "${r2out/.R./_R2.}" -Q -c -w 16 -D --trim_front1 6 --trim_front2 3 -h "${r1out/F.fq.gz/html}" -j "${r1out/F.fq.gz/json}"
done

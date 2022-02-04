#!/bin/sh

#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_fastp

cd $SLURM_SUBMIT_DIR

#specify conda environment
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate ddradQC 


#probably a better way to do this but whatever
#if data are not already demultiplexed then one can simply link to the fastq file containing all reads
dumbnumbers=("01" "02" "03" "04" "05" "06" "07" "08" "09" "DNT006" "PP56")

#for loop for the samples that do not conform to standard numbering scheme
for i in "${dumbnumbers[@]}"
do
	fastp -i /work/bs66/davidsonii_mapping/demuxed/PopF2_"$i".F.fq.gz -I /work/bs66/davidsonii_mapping/demuxed/PopF2_"$i".R.fq.gz -o PopF2_"$i".F.noadapters.fq.gz -O PopF2_"$i".R.noadapters.fq.gz -h output_fastp_"$i".html - j output_fastp_"$i".json -Q -L -w 10
done

#for loop for standard numbers.
for j in {10..84}
do
	fastp -i /work/bs66/davidsonii_mapping/demuxed/PopF2_"$j".F.fq.gz -I /work/bs66/davidsonii_mapping/demuxed/PopF2_"$j".R.fq.gz -o PopF2_"$j".F.noadapters.fq.gz -O PopF2_"$j".R.noadapters.fq.gz -h output_fastp_"$j".html - j output_fastp_"$j".json -Q -L -w 10
done

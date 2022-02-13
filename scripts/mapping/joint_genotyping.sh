#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -p wessinger-48core
#SBATCH --job-name=genotyping
#SBATCH --output=/work/bs66/davidsonii_mapping/mapping/log_outfiles/slurm-%j.out

cd $SLURM_SUBMIT_DIR

module load GATK/4.1.3.0
module load java

#important files and directories
genomefile="/work/bs66/davidsonii_mapping/davidsonii_genome.fasta"
vcffile="/work/bs66/davidsonii_mapping/mapping/genotyping/cohort_F2s.g.vcf.gz"
outdir="/work/bs66/davidsonii_mapping/mapping/genotyping"


gatk --java-options "-Xmx4g" GenotypeGVCFs -R $genomefile -V $vcffile -O $outdir/genotyped_cohort.vcf.gz

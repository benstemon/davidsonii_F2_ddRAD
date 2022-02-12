#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -p wessinger-48core
#SBATCH --job-name=genotyping
#SBATCH --output=/work/bs66/davidsonii_mapping/mapping/log_outfiles/slurm-%j.out

cd $SLURM_SUBMIT_DIR

module load GATK/4.1.3.0
module load java

#paths to important files and directories
genomefile="/work/bs66/davidsonii_mapping/davidsonii_genome.fasta"
infiledir="/work/bs66/davidsonii_mapping/mapping/haplotype_caller"
outfiledir="/work/bs66/davidsonii_mapping/mapping/genotyping"


gatk --java-options '-Xmx4g' CombineGVCFs\
 -R $genomefile\
 --variant $infiledir/PopF2_01.g.vcf.gz\
 --variant $infiledir/PopF2_02.g.vcf.gz\
 --variant $infiledir/PopF2_03.g.vcf.gz\
 --variant $infiledir/PopF2_04.g.vcf.gz\
 --variant $infiledir/PopF2_05.g.vcf.gz\
 --variant $infiledir/PopF2_06.g.vcf.gz\
 --variant $infiledir/PopF2_07.g.vcf.gz\
 --variant $infiledir/PopF2_08.g.vcf.gz\
 --variant $infiledir/PopF2_09.g.vcf.gz\
 --variant $infiledir/PopF2_10.g.vcf.gz\
 --variant $infiledir/PopF2_11.g.vcf.gz\
 --variant $infiledir/PopF2_12.g.vcf.gz\
 --variant $infiledir/PopF2_13.g.vcf.gz\
 --variant $infiledir/PopF2_14.g.vcf.gz\
 --variant $infiledir/PopF2_15.g.vcf.gz\
 --variant $infiledir/PopF2_16.g.vcf.gz\
 --variant $infiledir/PopF2_17.g.vcf.gz\
 --variant $infiledir/PopF2_18.g.vcf.gz\
 --variant $infiledir/PopF2_19.g.vcf.gz\
 --variant $infiledir/PopF2_20.g.vcf.gz\
 --variant $infiledir/PopF2_21.g.vcf.gz\
 --variant $infiledir/PopF2_22.g.vcf.gz\
 --variant $infiledir/PopF2_23.g.vcf.gz\
 --variant $infiledir/PopF2_24.g.vcf.gz\
 --variant $infiledir/PopF2_25.g.vcf.gz\
 --variant $infiledir/PopF2_26.g.vcf.gz\
 --variant $infiledir/PopF2_27.g.vcf.gz\
 --variant $infiledir/PopF2_28.g.vcf.gz\
 --variant $infiledir/PopF2_29.g.vcf.gz\
 --variant $infiledir/PopF2_30.g.vcf.gz\
 --variant $infiledir/PopF2_31.g.vcf.gz\
 --variant $infiledir/PopF2_32.g.vcf.gz\
 --variant $infiledir/PopF2_33.g.vcf.gz\
 --variant $infiledir/PopF2_34.g.vcf.gz\
 --variant $infiledir/PopF2_35.g.vcf.gz\
 --variant $infiledir/PopF2_36.g.vcf.gz\
 --variant $infiledir/PopF2_37.g.vcf.gz\
 --variant $infiledir/PopF2_38.g.vcf.gz\
 --variant $infiledir/PopF2_39.g.vcf.gz\
 --variant $infiledir/PopF2_40.g.vcf.gz\
 --variant $infiledir/PopF2_41.g.vcf.gz\
 --variant $infiledir/PopF2_42.g.vcf.gz\
 --variant $infiledir/PopF2_43.g.vcf.gz\
 --variant $infiledir/PopF2_44.g.vcf.gz\
 --variant $infiledir/PopF2_45.g.vcf.gz\
 --variant $infiledir/PopF2_46.g.vcf.gz\
 --variant $infiledir/PopF2_47.g.vcf.gz\
 --variant $infiledir/PopF2_48.g.vcf.gz\
 --variant $infiledir/PopF2_49.g.vcf.gz\
 --variant $infiledir/PopF2_50.g.vcf.gz\
 --variant $infiledir/PopF2_51.g.vcf.gz\
 --variant $infiledir/PopF2_52.g.vcf.gz\
 --variant $infiledir/PopF2_53.g.vcf.gz\
 --variant $infiledir/PopF2_54.g.vcf.gz\
 --variant $infiledir/PopF2_55.g.vcf.gz\
 --variant $infiledir/PopF2_56.g.vcf.gz\
 --variant $infiledir/PopF2_57.g.vcf.gz\
 --variant $infiledir/PopF2_58.g.vcf.gz\
 --variant $infiledir/PopF2_59.g.vcf.gz\
 --variant $infiledir/PopF2_60.g.vcf.gz\
 --variant $infiledir/PopF2_62.g.vcf.gz\
 --variant $infiledir/PopF2_63.g.vcf.gz\
 --variant $infiledir/PopF2_64.g.vcf.gz\
 --variant $infiledir/PopF2_65.g.vcf.gz\
 --variant $infiledir/PopF2_66.g.vcf.gz\
 --variant $infiledir/PopF2_67.g.vcf.gz\
 --variant $infiledir/PopF2_68.g.vcf.gz\
 --variant $infiledir/PopF2_69.g.vcf.gz\
 --variant $infiledir/PopF2_70.g.vcf.gz\
 --variant $infiledir/PopF2_71.g.vcf.gz\
 --variant $infiledir/PopF2_72.g.vcf.gz\
 --variant $infiledir/PopF2_73.g.vcf.gz\
 --variant $infiledir/PopF2_74.g.vcf.gz\
 --variant $infiledir/PopF2_75.g.vcf.gz\
 --variant $infiledir/PopF2_76.g.vcf.gz\
 --variant $infiledir/PopF2_77.g.vcf.gz\
 --variant $infiledir/PopF2_78.g.vcf.gz\
 --variant $infiledir/PopF2_79.g.vcf.gz\
 --variant $infiledir/PopF2_80.g.vcf.gz\
 --variant $infiledir/PopF2_81.g.vcf.gz\
 --variant $infiledir/PopF2_82.g.vcf.gz\
 --variant $infiledir/PopF2_83.g.vcf.gz\
 --variant $infiledir/PopF2_84.g.vcf.gz\
 --variant $infiledir/PopF2_DNT006.g.vcf.gz\
 --variant $infiledir/PopF2_PP56.g.vcf.gz\
 -O $outfiledir/cohort_F2s.g.vcf.gz


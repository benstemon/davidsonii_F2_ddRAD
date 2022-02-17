#these calculations are very quick and do not require jobs to be sumbitted
#can be performed easily on interactive node

module load vcftools/0.1.17

vcffile='/work/bs66/davidsonii_mapping/vcf_filtering/bestsnps_cohort_addheader.vcf'
outfile='/work/bs66/davidsonii_mapping/vcf_filtering'

#filter low GQ and DP genotypes 
vcftools --vcf $vcffile --minGQ 20 --minDP 4 --recode --recode-INFO-all --out $outfile/finalized_snps

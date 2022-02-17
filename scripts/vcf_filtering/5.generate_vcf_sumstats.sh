#these calculations are very quick and do not require jobs to be sumbitted
#can be performed easily on interactive node

module load vcftools/0.1.17

vcffile='/work/bs66/davidsonii_mapping/vcf_filtering/finalized_snps.recode.vcf'
outfile='/work/bs66/davidsonii_mapping/vcf_filtering/summary_outfiles/out_vcf'

#calculate allele frequency distributions
vcftools --vcf $vcffile --freq2 --out $outfile

#mean depth/individual
vcftools --vcf $vcffile --depth --out $outfile

#mean depth/site
vcftools --vcf $vcffile --site-mean-depth --out $outfile

#site quality
vcftools --vcf $vcffile --site-quality --out $outfile

#proportion missing data/individual
vcftools --vcf $vcffile --missing-indv --out $outfile

#proportion missing data/site
vcftools --vcf $vcffile --missing-site --out $outfile

#heterozygosity and inbreeding coefficient per individual
vcftools --vcf $vcffile --het --out $outfile


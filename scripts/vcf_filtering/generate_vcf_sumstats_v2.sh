#these calculations are very quick and do not require jobs to be sumbitted
#can be performed easily on interactive node

module load vcftools/0.1.17

infiledir='/work/bs66/davidsonii_mapping/mapping/vcf_filtering_no-phase/finalized_data_no-phase'
mkdir $infiledir/summary_outfiles


for i in $infiledir/*;
do
    #for moving throughout directories
    filehead="${i##*/}"
    
    #calculate allele frequency distributions
    vcftools --vcf $i --freq2 --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #mean depth/individual
    vcftools --vcf $i --depth --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #mean depth/site
    vcftools --vcf $i --site-mean-depth --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #site quality
    vcftools --vcf $i --site-quality --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #proportion missing data/individual
    vcftools --vcf $i --missing-indv --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #proportion missing data/site
    vcftools --vcf $i --missing-site --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #heterozygosity and inbreeding coefficient per individual
    vcftools --vcf $i --het --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #quality by depth 
    egrep -v "^#" $i | \
    cut -f 8 | \
    sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > $infiledir/summary_outfiles/out_vcf_$filehead.QD.txt
done
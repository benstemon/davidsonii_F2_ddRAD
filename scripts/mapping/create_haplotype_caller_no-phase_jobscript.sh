#full permanent paths for important files
genomefile="/work/bs66/davidsonii_mapping/davidsonii_genome.fasta"
headerfile="/work/bs66/davidsonii_mapping/mapping/haplotype_header.txt"

#full permanent paths for important directories
jobscriptdir="/work/bs66/davidsonii_mapping/mapping/log_outfiles"
infiledir="/work/bs66/davidsonii_mapping/mapping/RG"
outfiledir="/work/bs66/davidsonii_mapping/mapping/haplotype_caller_no-phase"


for i in $infiledir/*.bam;
do
    #naming variables
    readhead="${i##*/}"
    readhead="${readhead/.bam/}"
    
    #copy header file to jobscript
    cat $headerfile > $jobscriptdir/hapcaller_$readhead.sh
    
    #write HaplotypeCaller command to jobscript
    echo "gatk --java-options '-Xmx4g' HaplotypeCaller -R $genomefile -I $i -O $outfiledir/$readhead.g.vcf.gz -ERC GVCF --do-not-run-physical-phasing" >> $jobscriptdir/hapcaller_$readhead.sh
done


#full permanent paths for important files
genomefile="/work/bs66/davidsonii_mapping/davidsonii_genome.fasta"
headerfile="/work/bs66/davidsonii_mapping/mapping/mapping_header.txt"

#full permanent paths for important directories
infilepath="/work/bs66/davidsonii_mapping/preprocessing_v3/stacks_output"
jobscriptpath="/work/bs66/davidsonii_mapping/mapping/log_outfiles"
bamfilepath="/work/bs66/davidsonii_mapping/mapping/bamfiles"
rgpath="/work/bs66/davidsonii_mapping/mapping/RG"

#path to picard .jar file
picardpath="/share/apps/gcc/4.8.5/picard2018/picard.jar"


#main loop to create jobscripts
for r1in in $infilepath/*R1.1.fq.gz;
do
    #naming variables
    r2in="${r1in/R1.1/R2.2}"
    r1qz="${r1in##*/}"
    readhead="${r1qz/_R1.1.fq.gz/}"
    
    #copy header file to jobscript
    cat $headerfile > $jobscriptpath/mapping_$readhead.sh
    
    #write bwa command to jobscript (output = .sam file)
    echo "bwa mem $genomefile $r1in $r2in > $bamfilepath/$readhead.sam" >> $jobscriptpath/mapping_$readhead.sh
    
    #write samtools commands to jobscript (.sam -> unsorted.bam -> .bam)
    echo "samtools view -bS $bamfilepath/$readhead.sam > $bamfilepath/$readhead.unsorted.bam" >> $jobscriptpath/mapping_$readhead.sh
    echo "samtools sort $bamfilepath/$readhead.unsorted.bam > $bamfilepath/$readhead.bam" >> $jobscriptpath/mapping_$readhead.sh
    
    #write picard command to jobscript (Add or Replace Read Groups)
    echo "java -jar $picardpath AddOrReplaceReadGroups I=$bamfilepath/$readhead.bam O=$rgpath/$readhead.bam SO=coordinate RGID=SeqRUN# RGLB=$bamfilepath/$readhead.bam RGPL=illumina RGPU=$bamfilepath/$readhead.bam RGSM=$bamfilepath/$readhead.bam VALIDATION_STRINGENCY=LENIENT" >> $jobscriptpath/mapping_$readhead.sh
    
    #write final samtools command to jobscript (creates index file)
    echo "samtools index $rgpath/$readhead.bam $rgpath/$readhead.bai" >> $jobscriptpath/mapping_$readhead.sh
done


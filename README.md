# QC on Penstemon davidsonii F2 ddRAD data


## Working directly from pre-demuxed data
* data demultiplexed by Kate Ostevik using Sabre script (should link that here)


## Preprocessing: Approach 2 (fastp -> stacks, trimming 5' ends)


### Filter Illumina adapters, fix bases in read overlap, polyG trim, and trim restriction overhang with [fastp](https://github.com/OpenGene/fastp)

* Forward and Reverse reads both have low quality bases at the cut site. This is likely due to low base diversity. Because all reads have the same bases at the cut site, the sequencer has trouble reliably determining the bases.
* A few of the individuals DO have higher quality reads here (taken from different pools -- parents DNT006 and PP56). Here it is clear that the cut site is present, though there is an extra base ('C') on the 5' end of forward reads (the reverse reads are unaffected). I'm not sure why this happens, but it may be an artifact of the demultiplexing step, or maybe an issue with barcode ligation. To remedy this, trim the first 6 bp of forward reads for all reads (corresponding to restriction overhang for EcoRI 'AATTC' + the additional 'C' preceding this) and the first 3 bp of the reverse reads in fastp. All parameters used in fastp:
* See script [`run_fastp_v2.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/preprocessing/run_fastp_v2.sh)
* Illumina adapter trimming enabled by default
* PolyG tail trimming enabled by default
* Disable quality filtering
    - -Q
* Enable base correction in overlapped regions
    - -c (minimum 30 bp overlap = default)
* 16 threads
    - -w 16
* Do not perform deduplication; no need for duplication rate estimation
    - --dont_eval_duplication
* Trim first 6 bp of forward reads
    - --trim_front1 6
* Trim first 3 bp of reverse reads
    - --trim_front2 3

```shell
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
```

### Use process_radtags in stacks to filter low quality reads with sliding window approach

* Note that we will retain stacks in this pipeline because I like the sliding window approach to process_radtags better than the quality filtering options available in fastp.
* Requires three files:
    1. To generate job script headings, [`base_script.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/preprocessing/base_script.sh)
    2. To create job scripts, [`create_stacks_jobscripts.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/preprocessing/create_stacks_jobscripts.sh)
    3. To submit jobs, [`run_stacks_v2.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/preprocessing/run_stacks_v3.sh)

This pipeline uses stacks to do the following:
* Clean data, removing any read with an uncalled base
    - -c
* Discard reads with low quality scores
    - -q
* Sliding window of 15% of read length
    - -w 0.15 (default = 0.15)
* Read discarded if average score within sliding window drops below this value
    - -s 20 (default = 10)
* Drop reads less than 30 bp
    - --len_limit 30

Basic stacks syntax is as follows:
```
process_radtags --paired -1 $forward_read -2 $reverse_read -i gzfastq -o $out_dir -c -q -w 0.15 -s 20 --len_limit 30 --disable_rad_check
```

## Map with bwa
* Prior to mapping, the genome must be indexed. Perform indexing with:
```
bwa index davidsonii_genome.fasta
```
* The mapping process will require three separate files
    1. To generate job script headings: [`mapping_header.txt`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/mapping_header.txt)
    2. To create job scripts: [`create_mapping_jobscript.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/create_mapping_jobscript.sh)
    3. To submit jobs: [`masterscript_mapping.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/masterscript_mapping.sh)

### Jobscript creation and main pipeline
1. Local alignment in bwa
2. Convert .sam file from alignment to .bam file with samtools
3. Sort the .bam file with samtools
4. Assign reads to read groups with picard
5. Index final .bam file with samtools

```shell
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
```

## Call variants with GATK HaplotypeCaller and genotype with GenotypeGVCFs
HaplotypeCaller calls SNPs and indels simultaneously through local *de novo* assembly of haplotypes. It generates an intermediate GVCF which can then be used in GenotypeGVCFs (GVCF workflow) for sample genotyping. It works by defining active regions, determining haplotypes by assembling the active region, determining likelihoods of the haplotypes given read data, and then assigning sample genotypes using genotype likelihoods.


### Reference genome formatting
* Before calling variants, need to create a dictionary (reference.dict) file...
`gatk CreateSequenceDictionary -R davidsonii_genome.fasta`

* ... as well as an index file (reference.fai) for the reference genome:
`samtools faidx davidsonii_genome.fasta`

### Variant Calling
* need three files
1. haplotype_header.txt
2. create_haplotype_caller_jobscript.sh
3. masterscript_haplotype_caller.sh

`bash create_haplotype_caller_jobscript.sh`
`bash masterscript_haplotype_caller.sh`

gatk HaplotypeCaller -R genome.fasta -I input.bam -O output.g.vcf.gz --emitRefConfidence GVCF


these two can go together.
### Combine output into multi-sample GVCF with CombineGVCFs
### Joint Genotyping
gatk GenotypeGVCFs

### Variant Recalibration
gatk VQSR on the vcfs







# QC on Penstemon davidsonii F2 ddRAD data


## Working directly from pre-demuxed data
* data demultiplexed by Kate Ostevik using Sabre script (should link that here)


## Preprocessing: Approach 2 (fastp -> stacks, trimming 5' ends)
* Note: approach 1 (no trimming) failed, due to issues with enzyme cut sites (see below).


### Filter Illumina adapters, fix bases in read overlap, polyG trim, and trim restriction overhang with [fastp](https://github.com/OpenGene/fastp)
Forward and Reverse reads both have low quality bases at the cut site. This is likely due to low base diversity. Because all reads have the same bases at the cut site, the sequencer has trouble reliably determining the bases.<br />

A few of the individuals DO have higher quality reads here (taken from different pools -- parents DNT006 and PP56). Here it is clear that the cut site is present, though there is an extra base ('C') on the 5' end of forward reads (the reverse reads are unaffected). I'm not sure why this happens, but it may be an artifact of the demultiplexing step, or maybe an issue with barcode ligation. To remedy this, trim the first 6 bp of forward reads for all reads (corresponding to restriction overhang for EcoRI 'AATTC' + the additional 'C' preceding this) and the first 3 bp of the reverse reads in fastp.<br />

All parameters used in fastp:
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
    1. To generate job script headings, [`base_script.txt`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/preprocessing/base_script.txt)
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

```shell
gatk CreateSequenceDictionary -R davidsonii_genome.fasta
```

* ... as well as an index file (reference.fai) for the reference genome:

```shell
samtools faidx davidsonii_genome.fasta
```

### Variant Calling
* need three files
1. [`haplotype_header.txt`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/haplotype_header.txt)
2. [`create_haplotype_caller_jobscript.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/create_haplotype_caller_jobscript.sh)
3. [`masterscript_haplotype_caller.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/masterscript_haplotype_caller.sh)

Basic syntax of a jobscript is as follows:

```shell
gatk --java-options '-Xmx4g' HaplotypeCaller -R davidsonii_genome.fasta -I PopF2_01.bam -O PopF2_01.g.vcf.gz -ERC GVCF
```

### Combine output into multi-sample GVCF with CombineGVCFs
* See [`combine_gvcf.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/combine_gvcf.sh)
* Basic syntax: 
```
gatk --java-options '-Xmx4g' CombineGVCFs -R davidsonii_genome.fasta\
 --variant PopF2_01.g.vcf.gz\
 --variant PopF2_02.g.vcf.gz\
 --variant ...\
 -O cohort_F2s.g.vcf.gz
```

### Joint Genotyping with GenotypeGVCFs
* See [`joint_genotyping.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/joint_genotyping.sh)
```shell
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
```

## VCF filtering
### Filter for biallelic SNPs with MQ > 30 that represent fixed differences between the two parent species
* See [`1.filterVCF_v2.py`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/1.filterVCF_v2.py)
* Requirements:
    - .vcf.gz file created from previous step must be unzipped (gunzip)
    - Two parents, which are the last two individuals
* In-file parameters (requires editing script):
    - `nF2s` = number of F2s in .vcf
    - `minIndiv` = desired minimum number of individuals for which to filter a SNP
    - Specify infile and outfile parameters

```python
# -*- coding: utf-8 -*-

"""
Created on Wed Dec 20 18:04:58 2017

@author: carrie

this assumes the final two samples in the vcf are the parents.
"""

# modified by BWS 2/15/22
# Now disallows multi-bp variants in both ref and alt alleles
# Now disallows missing data in either parent
# Now includes parsing for phased genotypes

invcf = open('genotyped_cohort.vcf', 'rU')
outfile = open('filtered_cohort.vcf', 'w')

nF2s = 83
minIndiv = 50
nsamples = nF2s + 2
minMQ = 30

def skipHeader(line):
    cols = line.replace('\n','').split('\t')
    if len(cols) < 2:
        return 'header'
    elif cols[0] == '#CHROM':
        return 'header'
    else:
	return cols

def extractVCFfields(sampleData):
    """ Extract data from sample-level fields in VCF file """
    if sampleData != './.':
        fields = sampleData.split(':')
        if (len(fields) > 4):
            alleleDepths = fields[1].split(',')
            totDepth = fields[2]
            phreds = fields[4].split(',')
            return [totDepth, alleleDepths, phreds]
        else:
            return 'missing'
    else:
        return 'missing'

def findMapQual(infoField):
    info = infoField.split(';')
    listpos = len(info) - 1
    MQscore = 0

    while listpos > 0: # start searching for MQ from end of info field
        if (info[listpos].split('='))[0] == 'MQ': # found it!
            MQscore = float((info[listpos].split('='))[1])
            listpos = 0
        else:
            listpos -= 1 # keep searching
    return MQscore
###################################################################

for line in invcf:
    cols = skipHeader(line)
    if cols != 'header':
        MQscore = findMapQual(cols[7])
        #if MQ score beats threshold, collect reference and alternative call
        if MQscore >= minMQ:
            altBase = cols[4]
            refBase = cols[3]
            #if either the reference or alt call is > 1bp, do nothing
            if len(altBase) > 1 or len(refBase) > 1:
                pass
                
            else:
                calls = 0
                for j in range(9, 9+nsamples):
                    annot = extractVCFfields(cols[j])
                    if annot != 'missing':
                        calls += 1
                p1 = cols[nsamples-2 + 9]
                p2 = cols[nsamples-1 + 9]
                
                #updated code for missing and phased data:
                p1fields = p1.split(':')[0]
                p2fields = p2.split(':')[0]
                if any(x in p1fields for x in ('./.', '0/1', '0|1', p2fields)):
                    pass
                elif any(x in p2fields for x in ('./.', '0/1', '0|1')):
                    pass
                elif p1fields == '0|0' and p2fields == '0/0':
                    pass
                elif p1fields == '0/0' and p2fields == '0|0':
                    pass
                elif p1fields == '1|1' and p2fields == '1/1':
                    pass
                elif p1fields == '1/1' and p2fields == '1|1':
                    pass
                else:
                    if calls >= minIndiv:
                        outfile.write(line)

outfile.close()
```

### Find single best SNP per RADtag (most data + highest rare allele frequency)
* See [`2.find_best_snp_v2.py`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/2.find_best_snp_v2.py)
* This script just specifies which SNPs are selected. Directly after, use `3.thin_best_snp.py` (below) to extract SNPs to new .vcf file.
* Requirements:
    - filtered .vcf file (output from `1.filterVCF_v2.py`)
* In-file parameters (requires editing script):
    - `MinMinor` = Minimum number of individuals that can have minor allele for SNP to be selected
    - `plants` = Number of F2s + number of parents in .vcf
    - `radtaglength` = Minimum bp distance apart that selected SNPs can be (e.g., if data are 150bp SE data, this value should be 150)
    - Specify infile and outfile

```python
# This program thins to best snp per radtag
# written by JKK

# modified by Carrie
# assumes you've already filtered for quality

# modified by BWS 2/15/22
# Now includes parsing for phased genotypes
# Now includes additional parameter "radtaglength" -- minimum 

#parameters to modify
############################
MinMinor = 8 # min individuals that will have minor allele

vcf = open("filtered_cohort.vcf", "rU")
out2a = open('best_ids_cohort.txt', 'w')


plants = 83 + 2 #F2 + parents

radtaglength = 300#length of longest RadTags (e.g., 300 for 150bp PE reads)

############################


#do not modify:
last_scaff = ''
bestsnp = ''
lastpos = 0
cp = 0

bestlist = []

for line_idx, line in enumerate(vcf):
    cols = line.replace('\n', '').split('\t')
    scaff = cols[0]
    position = int(cols[1])
    ref_base = cols[3]
    alt_base = cols[4]
    if line_idx % 10000 == 0:
        print scaff, position

    mincc = 0

    if len(alt_base) == 1:
        datums = [0, 0, 0]
        for j in range(9, 9 + plants):
            if cols[j] != "./.":
                geno = cols[j].split(":")[0]
                if geno == "0/0" or geno == "0|0":
                    datums[0] += 1
                elif geno == "0/1" or geno == "0|1":
                    datums[1] += 1
                elif geno == "1/1" or geno == "1|1":
                    datums[2] += 1
                else:
                    print "strange genotype: ", geno
        mincc = min(datums[0] + datums[1], datums[2] + datums[1])
#this isn't quite finding the best SNP per radtag.
#It's finding the best SNP and ensuring there aren't others within 150 bp (or whatever).
#You could lose some data if cut sites are close to mapped SNPs but on different RADtags.
#The point though is to try to avoid obviously linked SNPs (eg same RADtag) and this does.
        if scaff != last_scaff or (position - lastpos) > radtaglength:
            if cp >= MinMinor:
                out2a.write(bestsnp)
            cp = mincc
            bestsnp = scaff + '\t' + cols[1] + '\n'
            last_scaff = scaff
            lastpos = position
        elif mincc > cp and position != lastpos:
            cp = mincc
            bestsnp = scaff + '\t' + cols[1] + '\n'

if cp >= MinMinor:
    out2a.write(bestsnp)  # last snp

out2a.close()
```

### Extract best SNP per RADtag
* See [`3.thin_best_snp.py`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/3.thin_best_snp.py)
* Requirements:
    - filtered .vcf file (output from `1.filterVCF_v2.py`)
    - best_ids .txt output from preceding step (`2.find_best_snp_v2.py`)
* In-file parameteres (requires editing script):
    - Specify best_ids .txt file, filtered .vcf file, and outfile

```python
# Code originally from Carrie Wessinger
# Modified by BWS 2/15/22

in2a = open('best_ids_cohort.txt', 'rU')#best_ids file
vcf = open("filtered_cohort.vcf", "rU")#filtered .vcf file
outfile = open("bestsnps_cohort.vcf", 'w')

bestlist = []
for line in in2a:
    cols = line.replace('\n', '').split('\t')
    bestlist.append([cols[0], cols[1]])

for line in vcf:
    cols = line.replace('\n', '').split('\t')
    tig = cols[0]
    pos = cols[1]
    for x in bestlist:
        if tig == x[0] and pos == x[1]:
            outfile.write(line)

outfile.close()
```






### Notes
* Something to consider for attempting parallelization.
```
--native-pair-hmm-threads 2 
```
* Could consider joining the GVCF merging step and the genotyping step.

* Consider Variant recalibration (gatk VQSR) on the joint gvcfs?


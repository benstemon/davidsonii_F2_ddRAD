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
        - Edit headings as necessary for computing cluster
    2. To create job scripts, [`create_stacks_jobscripts.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/preprocessing/create_stacks_jobscripts.sh)
        - Edit parameters in-text as necessary for your environment
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

**After editing parameters as needed, just: `bash create_stacks_jobscripts.sh` and `bash run_stacks_v2.sh` (submits jobscripts)**

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
        - Edit headings as necessary for computing cluster
    2. To create job scripts: [`create_mapping_jobscript.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/mapping/create_mapping_jobscript.sh)
        - Edit parameters as necessary for your environment
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

**After editing parameters as needed, just: `bash create_mapping_jobscript.sh` and `bash masterscript_mapping.sh` (submits jobscripts)**

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

**After editing necessary parameters in these files, just: `bash create_haplotpe_caller_jobscript.sh` and `bash masterscript_haplotype_caller.sh` (submits jobs)**

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
### Use vcftools to filter reads by genotype quality (GQ) and read depth (DP)
### NOTE: THIS STEP HAS BEEN REMOVED FROM THE PIPELINE!!!
* See [`filterby_gq_dp.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/filterby_gq_dp_v2.sh)
This script uses vcftools to do the following:
    - change genotypes with genotype quality (GQ) < 20 and/or filtered depth (DP) < 4 to missing (./.)




### Filter for biallelic SNPs with MQ > 30 that represent fixed differences between the two parent species
* See [`filterVCF_v4.py`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/filterVCF_v4.py)
* Requirements:
    - python 2.7
    - .vcf.gz file created from previous step must be unzipped (gunzip)
    - Two parents, which are the last two individuals
* In-file parameters (requires editing script):
    - `nF2s` = number of F2s in .vcf
    - `minIndiv` = desired minimum number of individuals for which to filter a SNP
        - I tested three missing data thresholds: minimums of 25, 42, and 50 individuals (70%, 50%, and 40% missing data/locus allowed, respectively)
    - Specify infile and outfile parameters


### Filter based on Hardy-Weinberg proportions
* [`See filterVCF.q.hwe_v2.py`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/filterVCF.q.hwe_v2.py)
* Only keep loci with HWE significance > 0.01
* Tested two minimum allele frequencies:
    - 0.4 ??? q ??? 0.6
    - 0.3 ??? q ??? 0.7
Note: after testing, only kept 0.3


### Find and extract single best SNP per RADtag (most data + highest rare allele frequency)
* See [`find_best_snp_v2.py`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/find_best_snp_v2.py) and [`thin_best_snp.py`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/thin_best_snp.py)
* `find_best_snp_v2.py` specifies which SNPs are selected. Directly after, use `thin_best_snp.py` to extract SNPs to new .vcf file.
* Requirements:
    - filtered .vcf file (output from `filterVCF_v3.py`)
* In-file parameters (requires editing script):
    - `MinMinor` = Minimum number of individuals that can have minor allele for SNP to be selected
        - I set this to 8 individuals
    - `plants` = Number of F2s + number of parents in .vcf
    - `radtaglength` = Minimum bp distance apart that selected SNPs can be
        - I set this to 300 bp
    - Specify infile and outfile parameters


### Use vcftools to calculate heterozygosity and inbreeding coefficients per individual, and other summary statistics
* Note that prior to this step I pasted the .vcf heading back onto the filtered output so vcftools would be able to recognize the file
* See [`generate_vcf_sumstats_v2.sh`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/generate_vcf_sumstats_v2.sh)


```shell
#these calculations are very quick and do not require jobs to be sumbitted
#can be performed easily on interactive node

module load vcftools/0.1.17

infiledir='/work/bs66/davidsonii_mapping/mapping/vcf_filtering/data'


for i in $infiledir/*;
do
    #for moving throughout directories
    filehead="${i##*/}"
    
    #calculate allele frequency distributions
    vcftools --vcf $i/finalized_snps_$filehead.vcf --freq2 --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #mean depth/individual
    vcftools --vcf $i/finalized_snps_$filehead.vcf --depth --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #mean depth/site
    vcftools --vcf $i/finalized_snps_$filehead.vcf --site-mean-depth --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #site quality
    vcftools --vcf $i/finalized_snps_$filehead.vcf --site-quality --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #proportion missing data/individual
    vcftools --vcf $i/finalized_snps_$filehead.vcf --missing-indv --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #proportion missing data/site
    vcftools --vcf $i/finalized_snps_$filehead.vcf --missing-site --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #heterozygosity and inbreeding coefficient per individual
    vcftools --vcf $i/finalized_snps_$filehead.vcf --het --out $infiledir/summary_outfiles/out_vcf_$filehead
    
    #quality by depth 
    egrep -v "^#" $i/finalized_snps_$filehead.vcf | \
    cut -f 8 | \
    sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > $infiledir/summary_outfiles/out_vcf_$filehead.QD.txt
done
```

### Calculate heterozygosity, read depth, etc. per individual
* See [`calc.sample.coverage.from.vcf_v3`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/scripts/vcf_filtering/calc.sample.coverage.from.vcf_v3.py)


### Visualize results of VCF filtering
* See [`davF2_filterVCF_results.pdf`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/results/davF2_filterVCF_results.pdf) for results and [`davF2_filterVCF_results.md`](https://github.com/benstemon/davidsonii_F2_ddRAD/blob/main/results/davF2_filterVCF_results.rmd) for code





## Create linkage map in Lep-Map3
### Convert .vcf to lepmap input
* See [`generate_lepmap_files_v2.py`]()
* Also requires plain .txt file with names of F2s in order -- see [`F2_namefile.txt`]()




### Notes
* Something to consider for attempting parallelization.
```
--native-pair-hmm-threads 2 
```
* Could consider joining the GVCF merging step and the genotyping step.

* Consider Variant recalibration (gatk VQSR) on the joint gvcfs?


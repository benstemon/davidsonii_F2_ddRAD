# QC on Penstemon davidsonii F2 ddRAD data


## Working directly from pre-demuxed data
* data demultiplexed by Kate Ostevik using Sabre script (should link that here)


## Preprocessing: Approach 2 (fastp -> stacks, trimming 5' ends)


### Filter Illumina adapters, fix bases in read overlap, polyG trim, deduplication, and trim restriction overhang with [fastp](https://github.com/OpenGene/fastp)

* Forward and Reverse reads both have low quality bases at the cut site. This is likely due to low base diversity. Because all reads have the same bases at the cut site, the sequencer has trouble reliably determining the bases.
* A few of the individuals DO have higher quality reads here (taken from different pools -- parents DNT006 and PP56). Here it is clear that the cut site is present, though there is an extra base ('C') on the 5' end of forward reads (the reverse reads are unaffected). I'm not sure why this happens, but it may be an artifact of the demultiplexing step, or maybe an issue with barcode ligation. To remedy this, trim the first 6 bp of forward reads for all reads (corresponding to restriction overhang for EcoRI 'AATTC' + the additional 'C' preceding this) and the first 3 bp of the reverse reads in fastp. All parameters used in fastp:
* Script available at `/davidsonii_F2_ddRAD/scripts/preprocessing/run_fastp_v2.sh'
* Illumina adapter trimming enabled by default
* PolyG tail trimming enabled by default
* Disable quality filtering
    - -Q
* Enable base correction in overlapped regions
    - -c (minimum 30 bp overlap = default)
* 16 threads
    - -w 16
* Evaluate duplication rate and remove duplicated reads
    - -D (default duplication calculation accuracy = 3)
* Trim first 6 bp of forward reads
    - --trim_front1 6
* Trim first 3 bp of reverse reads
    - --trim_front2 3

```
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
(Note that we will retain stacks in this pipeline because I like the sliding window approach to process_radtags better than the quality filtering options available in fastp.)

* Script available at `/davidsonii_F2_ddRAD/scripts/preprocessing/run_stacks_v2.sh'
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

```
#!/bin/sh

#SBATCH -N 1
#SBATCH -n 10 
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_stacks


cd $SLURM_SUBMIT_DIR


module load stacks/gcc/2.41


infilepath="/work/bs66/davidsonii_mapping/preprocessing_v2/fastp_outfiles"
outfilepath="/work/bs66/davidsonii_mapping/preprocessing_v2/stacks_output"


for r1in in $infilepath/*R1.fq.gz;
do
    r2in="${r1in/_R1./_R2.}"
    process_radtags --paired -1 $r1in -2 $r2in -i gzfastq -o $outfilepath -c -q -w 0.15 -s 20 --len_limit 30 --disable_rad_check
done
```

## Now moving on to alignment with BWA



# Testing different QC approaches to the davidsonii F2 ddRAD data

## Working directly from pre-demuxed data
* data demultiplexed by Kate Ostevik using Sabre script (should link that here)
* Two parent individuals (DNT and PP) appear to have proper RE overhang, with an additional 'C' base at the beginning. This pattern is not apparent in the F2 data (either raw or demux files)

### Preprocessing: Approach 1 (fastp -> stacks)
#### Filter Illumina adapters and fix bases in read overlap with [fastp](https://github.com/OpenGene/fastp)
* See `./davidsonii_F2_ddRAD/scripts/preprocessing/`
* Illumina adapter trimming enabled by default
* Disable quality filtering
    - -Q
* Remove reads of length <30 bp
    - --length_required 30
* Enable base correction in overlapped regions
    - -c (minimum 30 bp overlap = default)
* 16 threads
    - -w 16
* Evaluate duplication rate and remove duplicated reads
    - -D (default duplication calculation accuracy = 3)

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
    $fastpdir/./fastp -i "$r1in" -I "$r2in" -o "${r1out/PopF2_/PopF2_R1_}" -O "${r2out/PopF2_/PopF2_R2_}" -Q --length_required 30 -c -w 16 -D -h "${r1out/F.fq.gz/html}" -j "${r1out/F.fq.gz/json}"
done
```

#### Use process_radtags in stacks to rescue  restriction enzyme cutsite when possible and filter low quality reads with sliding window approach
* Clean data, removing any read with an uncalled base
    - -c
* Discard reads with low quality scores
    - -q
* Rescue RAD-Tags
    - -r
* Sliding window of 15% of read length
    - -w 0.15 (default = 0.15)
* Read discarded if average score within sliding window drops below this value
    - -s 10 (default = 10)
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


infilepath="/work/bs66/davidsonii_mapping/preprocessing_v1/fastp_outfiles"
outfilepath="/work/bs66/davidsonii_mapping/preprocessing_v1/stacks_output"


for r1in in $infilepath/PopF2_R1_*;
do
    r2in="${r1in/_R1_/_R2_}"
    process_radtags --paired -1 $r1in -2 $r2in -i gzfastq -o $outfilepath --renz_1 'ecoRI' --renz_2 'mspI' -c -q -r -w 0.15 -s 10 --len_limit 30
done
```

Interesting... there is a very high and very consistent amount of RAD cutsites not found (~80) for each sample.
Is there something I'm not understanding here?


### *ipyrad* testing... (includes strict filter for Illumina adapter)
#### referenced assembly, no RE recovery, no trimming:
#### referenced assembly, no RE recovery, trim first 5 bp of R1:
#### referenced assembly, no RE recovery, trim first 5 bp R1, first 3 bp R2:


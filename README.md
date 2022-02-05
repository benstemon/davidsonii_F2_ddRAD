# Testing different QC approaches to the davidsonii F2 ddRAD data

## Working directly from pre-demuxed data
* data demultiplexed by Kate Ostevik using Sabre script (should link that here)
* Two parent individuals (DNT and PP) appear to have proper RE overhang, with an additional 'C' base at the beginning. This pattern is not apparent in the F2 data (either raw or demux files)

### Preprocessing: Approach 1
#### Filter Illumina adapters and fix bases in read overlap with [fastp](https://github.com/OpenGene/fastp)
* See `./davidsonii_F2_ddRAD/scripts/preprocessing/`
* Illumina adapter trimming enabled by default
* low complexity filter enabled
    - -y (default 30% complexity)
* Remove reads of length <30 bp
    - --length_required 30
* Correct low-quality base calles in overlapping regions of paired-end reads
    - -c (minimum 30 bp overlap = default)
* 16 threads
    - -w 16
* Evaluate duplication rate and remove duplicated reads
    - -D (default duplication calculation accuracy = 3)
* Sliding window cutting enabled: 15% of read length, qs <15. This interferes with deduplication of data.
    - --cut_front
    - -W 22 (22 bp is ~15% of read length)
    - --cut_front_mean_quality 15

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
    $fastpdir/./fastp -i "$r1in" -I "$r2in" -o "${r1out/F.fq/F.postfastp.fq}" -O "${r2out/R.fq/R.postfastp.fq}" -y --length_required 30 -c -w 16 -D --cut_front -W 22 --cut_front_mean_quality 15 -h "${r1out/F.fq.gz/html}" -j "${r1out/F.fq.gz/json}"
done
```

## Next step is to... run this through bwa I suppose?
## Need some kind of metric for seeing if this worked ok. So... ?
## investigate these html files I suppose...


### *ipyrad* testing... (includes strict filter for Illumina adapter)
#### referenced assembly, no RE recovery, no trimming:
#### referenced assembly, no RE recovery, trim first 5 bp of R1:
#### referenced assembly, no RE recovery, trim first 5 bp R1, first 3 bp R2:








# this didn't work great, so probably not going to do it (didn't like file type, and I don't know enough about it..."

#### Fix the restriction site using python script 'recoverRE.py'
* (original code from [Matt Gibson](https://github.com/gibsonMatt/pimpGEA/blob/master/scripts/recoverREsite/recoverRE.py)
* See `./davidsonii_F2_ddRAD/scripts/recover_RE/`

```
python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out 



#Read 1
python recoverRE.py --truesite TGCAG --startpos 6 --endpos 10 --out /N/dc2/projects/gibsonTomato/pimpGEA/data/recoverREsite/poolA/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R1_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/pimpGEA/data/adapter_clean_rawdata/poolA/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R1_001_noadapters.fastq

#Read 2
python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out /N/dc2/projects/gibsonTomato/pimpGEA/data/recoverREsite/poolA/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R2_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/pimpGEA/data/adapter_clean_rawdata/poolA/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R2_001_noadapters.fastq

#...done for all index libraries
```











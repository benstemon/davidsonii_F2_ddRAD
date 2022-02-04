# Testing different QC approaches to the davidsonii F2 ddRAD data

## Working directly from pre-demuxed data
* data demultiplexed by Kate Ostevik using Sabre script (should link that here)
* Two parent individuals (DNT and PP) appear to have proper RE overhang, with an additional 'C' base at the beginning. This pattern is not apparent in the F2 data (either raw or demux files)

### Approach 1: fastp + fix restriction enzyme site

#### Filter Illumina adapters and fix bases in read overlap with [fastp](https://github.com/OpenGene/fastp)
* See `./davidsonii_F2_ddRAD/scripts/filter_adapters/`

```
#!/bin/sh

#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_fastp

cd $SLURM_SUBMIT_DIR

#specify conda environment
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate ddradQC 


#probably a better way to do this but whatever
#if data are not already demultiplexed then one can simply link to the fastq file containing all reads
dumbnumbers=("01" "02" "03" "04" "05" "06" "07" "08" "09" "DNT006" "PP56")

#for loop for the samples that do not conform to standard numbering scheme
for i in "${dumbnumbers[@]}"
do
	fastp -i /work/bs66/davidsonii_mapping/demuxed/PopF2_"$i".F.fq.gz -I /work/bs66/davidsonii_mapping/demuxed/PopF2_"$i".R.fq.gz -o PopF2_"$i".F.noadapters.fq.gz -O PopF2_"$i".R.noadapters.fq.gz -h output_fastp_"$i".html - j output_fastp_"$i".json -Q -L -w 10
done

#for loop for standard numbers.
for j in {10..84}
do
	fastp -i /work/bs66/davidsonii_mapping/demuxed/PopF2_"$j".F.fq.gz -I /work/bs66/davidsonii_mapping/demuxed/PopF2_"$j".R.fq.gz -o PopF2_"$j".F.noadapters.fq.gz -O PopF2_"$j".R.noadapters.fq.gz -h output_fastp_"$j".html - j output_fastp_"$j".json -Q -L -w 10
done
```

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




### *ipyrad* testing... (includes strict filter for Illumina adapter)
#### referenced assembly, no RE recovery, no trimming:
#### referenced assembly, no RE recovery, trim first 5 bp of R1:
#### referenced assembly, no RE recovery, trim first 5 bp R1, first 3 bp R2:




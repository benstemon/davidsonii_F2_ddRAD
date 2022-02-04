# Penstemon davidsonii processing F2 ddRAD data

## Pre-process data
* Reads were already demultiplexed with [sabre](https://github.com/najoshi/sabre). This removes barcodes (but presumably still leaves adapter cut sites?)
* This lead to the identification of an issue with the data. There appears to be a large degree of misspecification with the ligation of adapters, such that they have ligated to anything in the digest, rather than specific to restriction enzyme cut-sites. There is likely also a problem with sequencing error due to low base diversity.
* It's possible to use FASTP to (1) change base pairs to avoid mapping errors [(Gibson and Moyle, 2020)](https://onlinelibrary.wiley.com/doi/pdf/10.1111/mec.15477)

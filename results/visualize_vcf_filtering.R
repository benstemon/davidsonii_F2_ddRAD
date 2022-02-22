#using tidyverse for manipulating output
library(tidyverse)

#Taking some pointers from Mark Ravinet & Joana Meier and materials from their Physalia course:
#https://speciationgenomics.github.io/filtering_vcfs/



###
setwd('~/workdir/data/summary_outfiles')

#going to use the following filtered data set:
#min 50; hwe 0.3 ≤ q ≤ 0.7; RADtaglength 300bp

prefix <- 'out_vcf_min50_hwe.30_300bp'
###

#The current data set has 1427 SNPs. This is after stringent quality filtering.
# minimum Genotype quality is 20 (99% accuracy)
# minimum allele depth is 4
# minimum Mapping Quality is 30
# no more than ~40% missing data (50 or more individuals must be present)
# allele frequencies must be in hardy-weinberg proportions
# Single SNP per 300 bp
# at least 8 individuals with minor allele



#
##QUALITY BY DEPTH
#
#plotting QD score (qual by depth) (not extracted from vcftools)
#recommendation is to filter out QD < 2
t <- read_table(paste(prefix, '.QD.txt', sep=''), col_names = 'QD')
ggplot(t, aes(QD)) +
  geom_density(fill = "orange", colour = "black", alpha = 0.3) +
  xlim(0,40) +
  ggtitle("Quality by Depth") +
  theme_classic()

#looks very good. We have no low quality sites
#after filtering, nothing with QD<5
length(which(t<5))



#
##DEPTH OF COVERAGE
#
#Higher coverage is better, obviously
#but, reads with too high coverage could be mapping/assembly errors and/or repetitive regions
#Ravinet & Meier suggest good "rule of thumb" is filtering max depth > 2x mean depth
#But I have seen less stringent filters elsewhere
t <- read_delim(paste(prefix, '.ldepth.mean', sep = ''), delim = '\t')
ggplot(t, aes(MEAN_DEPTH)) + 
  geom_density(fill = "blue", colour = "black", alpha = 0.3) +
  xlim(0,150) +
  ggtitle("Quality by Depth") +
  theme_classic()

#This looks pretty good.
#if we look for the proportion of reads > 2x mean depth...
length(which(t$MEAN_DEPTH > mean(t$MEAN_DEPTH)*2))/nrow(t)
#12.3% are higher than 2x mean. But none are particularly high coverage.
#given this is ddrad data, nothing here screams mapping error to me.

#and nothing with particularly low coverage:
length(which(t$MEAN_DEPTH < 10))/nrow(t)
length(which(t$MEAN_DEPTH < 5))/nrow(t)


#
##MISSING DATA
#
t <- read_delim(paste(prefix, '.lmiss', sep = ''), delim = '\t')
ggplot(t, aes(F_MISS)) + 
  geom_density(fill = "red", colour = "black", alpha = 0.3) +
  xlim(0,1)
#looks how we would expect:
#we filtered for no more than 33 individuals with missing data (~40%)


#
##MINOR ALLELE FREQUENCY
#
t <- read_delim(paste(prefix, '.frq', sep = ''), delim = '\t',
                col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
#this is just allele freqs. Need to calc MAF.
t$maf <- t %>% select(a1, a2) %>% apply(1, function(z) min(z))
ggplot(t, aes(maf)) + 
  geom_density(fill = "black", colour = "black", alpha = 0.3) +
  xlim(0,0.5)
#Again, we filtered this so that minor allele frequency is always > 0.3


#
##HETEROZYGOSITY
#
library(gridExtra)

cov <- read.delim('~/workdir/vcf.sample.coverage.txt', header = T)
F2s <- droplevels(filter(cov, sample != 'PopF2_DNT006' & sample != 'PopF2_PP56'))

#Total sites
ggplot(F2s, aes(sites)) +
  geom_density(fill = "grey", col = "black") +
  xlim(500,1500) +
  ggtitle("Total sites") +
  theme_classic()

#Median genotypes per site
ggplot(F2s, aes(median_genos_persite)) +
  geom_density(fill = "grey", col = "black") +
  xlim(70,85) +
  ggtitle("Median genotypes per site") +
  theme_classic()

#Median depth per site
ggplot(F2s, aes(median_depth_persite)) +
  geom_density(fill = "grey", col = "black") +
  xlim(0,60) +
  ggtitle("Median depth per site") +
  theme_classic()

#Heterozygosity by number of sites
ggplot(F2s, aes(x = sites, y = het_freq)) +
  geom_point() +
  theme_classic() +
  ggtitle("Heterozygosity by number of sites")



#Notes:
#My depth/site is much higher than Carrie's
#This is presumably because I filtered GQ < 20
#With a less stringent filter there we could return much more data







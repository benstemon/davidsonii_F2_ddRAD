# -*- coding: utf-8 -*-

"""
Created on Wed Dec 20 18:04:58 2017

@author: carrie

this assumes the final two samples in the vcf are the parents.
"""

# modified by BWS
# since v2:
# Now disallows multi-bp variants in both ref and alt alleles
# Now disallows missing data in either parent (problem may be unique to HaplotypeCaller?)
# Now includes parsing for phased genotypes
# since v3: includes parsing for missing variants in parents on phased haplotype

invcf = open('/work/bs66/davidsonii_mapping/mapping/vcf_filtering/filtered_genotyped_cohort.vcf.recode.vcf', 'rU')
outfile = open('/work/bs66/davidsonii_mapping/mapping/vcf_filtering/biallelic_filteredMQ.vcf', 'w')

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
                if any(x in p1fields for x in ('./.', '.|.', '0/1', '0|1', p2fields)):
                    pass
                elif any(x in p2fields for x in ('./.', '.|.', '0/1', '0|1')):
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

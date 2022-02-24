# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 15:04:53 2017

@author: carrie

getting basic info from vcf files

Assumes a filtered vcf!
"""
#updated by BWS 2/21/22
#now compatible with phased data from HaplotypeCaller output

from numpy import median
from numpy import mean
from numpy import sum
from scipy.stats import chi2

in1 = open('biallelic_fixed-parents.vcf', 'rU')
outfile = open('biallelic_fixed-parents_hwe.30.vcf', 'w')

nF2s = 83
nsamples = nF2s + 2 
minq = 0.3
maxq = 0.7
minHWE = 0.01

for line in in1:
    cols = line.replace('\n', '').split('\t')

    # Pass over metadata and header lines
    if len(cols) < 2:
        pass
    elif cols[0] == "#CHROM":
        pass

    else:
        tig = cols[0]
        pos = cols[1]
        totREF = 0
        totALT = 0
        calls = 0
        reads = []
        RRcount = 0
        RAcount = 0
        AAcount = 0

        for j in range(9, 9+nsamples):
        #update here for HaplotypeCaller output
            if cols[j].split(':')[0] != './.' and cols[j].split(':')[0] != '.|.':
                calls += 1
                fields = cols[j].split(':')  # e.g., 0/0:67,1:67:99:0,160,2179

                if fields[0] != './.' and fields[1] != './.' and fields[2] != './.':
                    alleleDepths = fields[1].split(',')
                    refDepth = int(alleleDepths[0])
                    altDepth = int(alleleDepths[1])
                    reads.append(float(refDepth + altDepth))

#another update from HaplotypeCaller
#not all fields are same length (if haplotype is phased there are more columns)
                    if len(fields) == 5:
                        phreds = fields[4].split(',')
                    elif len(fields) == 8:
                        phreds = fields[6].split(',')
                    if phreds[0] == '0':
                        RRcount += 1
                    #changed statement so only counts hets if it's the only 0 phred score
                    elif phreds[1] == '0' and phreds[0]!= '0' and phreds[2]!= '0':
                        RAcount += 1
                    elif phreds[2] == '0':
                        AAcount += 1

        meanRPI = mean(reads)
        sumRPI = sum(reads)
        freqRR = float(RRcount) / calls
        freqRA = float(RAcount) / calls
        freqAA = float(AAcount) / calls
        freqR = freqRR + 0.5 * freqRA
        freqA = freqAA + 0.5 * freqRA
        expRR = calls * freqR * freqR
        expRA = calls * 2 * freqA * freqR
        expAA = calls * freqA * freqA

        if expRR == 0:
            termRR = 0
        else:
            termRR = (((RRcount - expRR) ** 2) / expRR)

        if expRA == 0:
            termRA = 0
        else:
            termRA = (((RAcount - expRA) ** 2) / expRA)

        if expAA == 0:
            termAA = 0
        else:
            termAA = (((AAcount - expAA) ** 2) / expAA)

        chisq = termRR + termRA + termAA

        pvalue = 1.0 - chi2.cdf(chisq, 1)

        if pvalue >= minHWE and (minq <= freqR <= maxq):
            outfile.write(line)

        

outfile.close()


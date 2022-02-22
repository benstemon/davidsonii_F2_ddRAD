'''
author: Carrie Wessinger

this finds site coverage per sample

input is vcf file and list of samples
'''
#Modified by BWS 2/22/2022
#extractVCFfields function now works for phased + unphased data types


import numpy as np


vcf = open('finalized_snps_min50_hwe.30_300bp.vcf', 'rU')
outfile = open('vcf.sample.coverage.txt', 'w')

#this is a list of sample names (1/row)
labelfile = open('namefile.txt', 'rU')

nindivs = 83 + 2



##################################
def skipHeader(line):
    cols = line.replace('\n','').split('\t')
    if len(cols) < 2:
        return 'header'
    elif cols[0] == '#CHROM':
        return 'header'
    else:
        return cols

#updated this function for variable length columns (phased vs. not)
def extractVCFfields(sampleData):
    """ Extract data from sample-level fields in VCF file """
    if any(y in sampleData.split(':')[0] for y in ('./.', '.|.')):
        return 'missing'
    else:
        fields = sampleData.split(':')
        if (len(fields) == 5):
            alleleDepths = fields[1].split(',')
            totDepth = fields[2]
            phreds = fields[4].split(',')
            return [totDepth, alleleDepths, phreds]
        elif (len(fields) == 8):
            alleleDepths = fields[1].split(',')
            totDepth = fields[2]
            phreds = fields[6].split(',')
            return [totDepth, alleleDepths, phreds]
        else:
            return 'missing'

##################################


# read in labels for samples
# store in list = "labs"

labs = []
for line in labelfile:
    cols = line.replace('\n', '').split('\t')
    labs.append([cols[0], 0, [], [], 0])  # ['name', nsites, [genos_per_site], [depth_per_site], nhetsites]

lines = 0
for line in vcf:
    cols = skipHeader(line)
    if cols != 'header':
        lines += 1
        ngenos = 0
        for j in range(9, nindivs+9):
            if cols[j].split(':')[0] != './.':
                ngenos += 1
        for j in range(9, nindivs+9):
            annot = extractVCFfields(cols[j])
            if annot != 'missing':
                depth = annot[0]
                labs[j-9][1] += 1
                labs[j-9][2].append(ngenos)
                labs[j-9][3].append(float(depth))
                if annot[2][1] == '0':
                    labs[j-9][4] += 1
                else:
                    pass
            else:
                pass


#write column names
outfile.write('sample\tsites\tmedian_genos_persite\tmedian_depth_persite\thet_freq\n')


for l in labs:
    outfile.write(l[0] + '\t' + str(l[1]) + '\t' + str(np.median(l[2])) + '\t' + str(np.median(l[3])) + '\t' + str(float(l[4])/float(l[1])) + '\n')
outfile.close()

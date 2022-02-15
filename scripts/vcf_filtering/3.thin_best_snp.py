# Code originally from Carrie Wessinger
# Modified by BWS 2/15/22

in2a = open('best_ids_cohort.txt', 'rU')#best_ids file
vcf = open("filtered_cohort.vcf", "rU")#filtered .vcf file
outfile = open("bestsnps_cohort.vcf", 'w')

bestlist = []
for line in in2a:
    cols = line.replace('\n', '').split('\t')
    bestlist.append([cols[0], cols[1]])

for line in vcf:
    cols = line.replace('\n', '').split('\t')
    tig = cols[0]
    pos = cols[1]
    for x in bestlist:
        if tig == x[0] and pos == x[1]:
            outfile.write(line)

outfile.close()


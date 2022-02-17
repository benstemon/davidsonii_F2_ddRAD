# This program thins to best snp per radtag
# written by JKK

# modified by Carrie
# assumes you've already filtered for quality

# modified by BWS 2/15/22
# Now includes parsing for phased genotypes
# Now includes additional parameter "radtaglength" -- minimum 

#parameters to modify
############################
MinMinor = 8 # min individuals that will have minor allele

vcf = open("/work/bs66/davidsonii_mapping/mapping/vcf_filtering/biallelic_filteredMQ.vcf", "rU")
out2a = open('/work/bs66/davidsonii_mapping/mapping/vcf_filtering/best_ids_cohort.txt', 'w')


plants = 83 + 2 #F2 + parents

radtaglength = 300#length of longest RadTags (e.g., 300 for 150bp PE reads)

############################


#do not modify:
last_scaff = ''
bestsnp = ''
lastpos = 0
cp = 0

bestlist = []

for line_idx, line in enumerate(vcf):
    cols = line.replace('\n', '').split('\t')
    scaff = cols[0]
    position = int(cols[1])
    ref_base = cols[3]
    alt_base = cols[4]
    if line_idx % 10000 == 0:
        print scaff, position

    mincc = 0

    if len(alt_base) == 1:
        datums = [0, 0, 0]
        for j in range(9, 9 + plants):
            if cols[j] != "./." or ".|.":
                geno = cols[j].split(":")[0]
                if geno == "0/0" or geno == "0|0":
                    datums[0] += 1
                elif geno == "0/1" or geno == "0|1":
                    datums[1] += 1
                elif geno == "1/1" or geno == "1|1":
                    datums[2] += 1
                else:
                    print "strange genotype: ", geno
        mincc = min(datums[0] + datums[1], datums[2] + datums[1])
#this isn't quite finding the best SNP per radtag.
#It's finding the best SNP and ensuring there aren't others within 150 bp (or whatever).
#You could lose some data if cut sites are close to mapped SNPs but on different RADtags.
#The point though is to try to avoid obviously linked SNPs (eg same RADtag) and this does.
        if scaff != last_scaff or (position - lastpos) > radtaglength:
            if cp >= MinMinor:
                out2a.write(bestsnp)
            cp = mincc
            bestsnp = scaff + '\t' + cols[1] + '\n'
            last_scaff = scaff
            lastpos = position
        elif mincc > cp and position != lastpos:
            cp = mincc
            bestsnp = scaff + '\t' + cols[1] + '\n'

if cp >= MinMinor:
    out2a.write(bestsnp)  # last snp

out2a.close()



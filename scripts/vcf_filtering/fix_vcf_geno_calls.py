import csv
import sys


#infile is tab-delimited vcf file (no header)
#first column = marker name
#subsequent columns = individuals

infile = sys.argv[1]
outfile = sys.argv[2]

#
data_dict = {}


#convert vcf calls into ASmap format
with open(infile, 'r') as indata:
    for linecount, line in enumerate(indata):
        snplist = []
        marker = line.split('\t')
        for column in marker[9:len(marker)]:
            geno = column.split(':')
            if geno[0] == './.':
                outcall = './.'
            else:
                phreds = geno[4].split(',')
                if phreds[0] == '0' and phreds[1] == '0':
                    outcall = './.'
                elif phreds[1] == '0' and phreds[2] == '0':
                    outcall = './.'
                elif phreds[0] == '0' and phreds[2] == '0':
                    outcall = './.'
                else:
                    outcall = geno
            snplist.append(outcall)
        data_dict[linecount] = marker[0:9]+snplist


#write output file
with open(outfile, "wb") as f:
    csv_output = csv.writer(f)
    for key in data_dict.keys():
        csv_output.writerow([key] + data_dict[key])




egrep -v "^#" finalized_snps.vcf | \
cut -f 8 | \
sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > summary_outfiles/out_vcf.QD.txt

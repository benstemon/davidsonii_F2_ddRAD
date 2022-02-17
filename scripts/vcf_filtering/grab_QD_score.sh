egrep -v "^#" finalized_snps.recode.vcf | \
cut -f 8 | \
sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > summary_outfiles/QD.txt

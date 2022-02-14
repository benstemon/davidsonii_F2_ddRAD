infilepath="/work/bs66/davidsonii_mapping/preprocessing_v3/fastp_outfiles"
outfilepath="/work/bs66/davidsonii_mapping/preprocessing_v3/stacks_output"
jobscriptpath="/work/bs66/davidsonii_mapping/preprocessing_v3/stacks_jobscripts"

for r1in in $infilepath/*R1.fq.gz;
do
    r1name="${r1in##*/}"
    cat base_script.txt > $jobscriptpath/$r1name
    r2in="${r1in/R1/R2}"
    echo "process_radtags --paired -1 $r1in -2 $r2in -i gzfastq -o $outfilepath -c -q -w 0.15 -s 20 --len_limit 30 --disable_rad_check" >> $jobscriptpath/$r1name
    mv $jobscriptpath/$r1name $jobscriptpath/stacks_"${r1name/_R1.fq.gz/.sh}"
done

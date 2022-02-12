#script to batch submit jobs
jobscriptpath="/work/bs66/davidsonii_mapping/mapping/log_outfiles"

for i in $jobscriptpath/mapping_*.sh;
do
    sbatch $i
done

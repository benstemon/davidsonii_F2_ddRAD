#script to batch submit jobs
jobscriptpath="/work/bs66/davidsonii_mapping/mapping/log_outfiles"

for i in $jobscriptpath/*.sh;
do
    sbatch $i
done

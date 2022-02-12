#script to batch submit jobs for HaplotypeCaller
jobscriptpath="/work/bs66/davidsonii_mapping/mapping/log_outfiles"

for i in $jobscriptpath/hapcaller_*.sh;
do
    sbatch $i
done

#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0635
#PBS -N tpstats_lpanamensis-hpgl0635 -l mem=1gb -l walltime=00:10:00 -l ncpus=1
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0635/outputs/qsub/tpstats_lpanamensis-hpgl0635.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0635/scripts/tpstats_lpanamensis-hpgl0635.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0635 || exit


## This is a stupidly simple job to collect tophat alignment statistics.


if [ ! -r outputs/tophat_stats.csv ]; then
  echo "basename,species,original_reads,aligned_reads,failed_reads,rpm,count_table" > outputs/tophat_stats.csv
fi
bamtools stats < outputs/tophat_lpanamensis/accepted_hits.bam 2>outputs/tophat_lpanamensis/accepted_hits.bam.stats 1>&2 && bamtools stats < outputs/tophat_lpanamensis/unmapped.bam 2>outputs/tophat_lpanamensis/unmapped.bam.stats 1>&2
original_reads_tmp=$(grep "^Input Reads" outputs/hpgl0635-trimmed-trimomatic.out | awk '{print $3}' | sed 's/ //g')
original_reads=${original_reads_tmp:-0}
reads_tmp=$(grep "^reads_in " outputs/tophat_lpanamensis/prep_reads.info | awk -F= '{print $2}' | sed 's/ //g')
reads=${reads_tmp:-0}
aligned_tmp=$(grep "^Total reads" outputs/tophat_lpanamensis/accepted_hits.bam.stats | awk '{print $3}' | sed 's/ .*//g')
aligned=${aligned_tmp:-0}
failed_tmp=$(grep "^Total reads" outputs/tophat_lpanamensis/unmapped.bam.stats | awk '{print $3}' | sed 's/ .*//g')
failed=${failed_tmp:-0}
rpm_tmp=$(perl -e "printf(1000000 / ${aligned})" 2>/dev/null)
rpm=${rpm_tmp:-0}
stat_string=$(printf "hpgl0635-trimmed,lpanamensis,%s,%s,%s,%s,%s,accepted_hits.count.xz" "${original_reads}" "${reads}" "${aligned}" "${failed}" "$rpm")
echo "$stat_string" >> outputs/tophat_stats.csv

## The following lines give status codes and some logging
echo $? > outputs/status/tpstats_lpanamensis-hpgl0635.status
echo "###Finished ${PBS_JODID} tpstats_lpanamensis-hpgl0635.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


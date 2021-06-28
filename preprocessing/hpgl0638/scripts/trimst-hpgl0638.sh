#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0638
#PBS -N trimst-hpgl0638 -l mem=1gb -l walltime=00:10:00 -l ncpus=1
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0638/outputs/qsub/trimst-hpgl0638.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0638/scripts/trimst-hpgl0638.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0638 || exit

## This is a stupidly simple job to collect trimomatic statistics

if [ ! -r outputs/trimomatic_stats.csv ]; then
  echo "total_reads,surviving_both,surviving_forward,surviving_reverse,dropped_reads" > outputs/trimomatic_stats.csv
fi
total_reads_tmp=$(grep "^Input Read Pairs" outputs/hpgl0638-trimomatic.out | awk '{print $4}')
total_reads=${total_reads_tmp:-0}
surviving_both_tmp=$(grep "^Input Read Pairs" outputs/hpgl0638-trimomatic.out | awk '{print $7}')
surviving_both=${surviving_both_tmp:-0}
surviving_forward_tmp=$(grep "^Input Read Pairs" outputs/hpgl0638-trimomatic.out | awk '{print $12}')
surviving_forward=${surviving_forward_tmp:-0}
surviving_reverse_tmp=$(grep "^Input Read Pairs" outputs/hpgl0638-trimomatic.out | awk '{print $17}')
surviving_reverse=${surviving_reverse_tmp:-0}
dropped_reads_tmp=$(grep "^Input Read Pairs" outputs/hpgl0638-trimomatic.out | awk '{print $20}')
dropped_reads=${dropped_reads_tmp:-0}

stat_string=$(printf "hpgl0638,%s,%s,%s,%s,%s" "${total_reads}" "${surviving_both}" "${surviving_forward}" "${surviving_reverse}" "${dropped_reads}")
echo "$stat_string" >> outputs/trimomatic_stats.csv

## The following lines give status codes and some logging
echo $? > outputs/status/trimst-hpgl0638.status
echo "###Finished ${PBS_JODID} trimst-hpgl0638.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


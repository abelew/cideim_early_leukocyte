#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/lab/nelsayed/raw_data/ade/HPGL0241/processed
#PBS -N unxz-HPGL0241 -l mem=4gb -l walltime=18:00:00 -l ncpus=1
#PBS -o /cbcb/lab/nelsayed/raw_data/ade/HPGL0241/processed/outputs/qsub/unxz-HPGL0241.qsubout -j oe -V -m n
echo "####Started /cbcb/lab/nelsayed/raw_data/ade/HPGL0241/processed/scripts/unxz-HPGL0241.sh at $(date)" >> outputs/log.txt
cd /cbcb/lab/nelsayed/raw_data/ade/HPGL0241/processed || exit

xz -f -d HPGL0241_R1_filtered_adaptertrim.fastq.xz HPGL0241_R2_filtered_adaptertrim.fastq.xz

## The following lines give status codes and some logging
echo $? > outputs/status/unxz-HPGL0241.status
echo "###Finished ${PBS_JODID} unxz-HPGL0241.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


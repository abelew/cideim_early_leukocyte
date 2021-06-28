#PBS -V -S /usr/bin/bash -q workstation
#PBS -d /cbcb/lab/nelsayed/raw_data/ade/HPGL0242/processed
#PBS -N fqc-hpgl0242 -l mem=6gb -l walltime=10:00:00 -l ncpus=8
#PBS -o /cbcb/lab/nelsayed/raw_data/ade/HPGL0242/processed/outputs/qsub/fqc-hpgl0242.qsubout -j oe -V -m n
echo "####Started /cbcb/lab/nelsayed/raw_data/ade/HPGL0242/processed/scripts/fqc-hpgl0242.sh at $(date)" >> outputs/log.txt
cd /cbcb/lab/nelsayed/raw_data/ade/HPGL0242/processed || exit

## This FastQC run is against rnaseq data and is used for
## an initial estimation of the overall sequencing quality.
mkdir -p outputs/fastqc && \
  fastqc --extract -o outputs/fastqc hpgl0242_forward.fastq.gz hpgl0242_reverse.fastq.gz \
  2>outputs/fastqc.out 1>&2

## The following lines give status codes and some logging
echo $? > outputs/status/fqc-hpgl0242.status
echo "###Finished ${PBS_JODID} fqc-hpgl0242.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


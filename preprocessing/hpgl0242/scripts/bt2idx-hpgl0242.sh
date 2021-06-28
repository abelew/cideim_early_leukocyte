#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/lab/nelsayed/raw_data/ade/HPGL0242/processed
#PBS -N bt2idx-hpgl0242 -l mem=6gb -l walltime=10:00:00 -l ncpus=4
#PBS -o /cbcb/lab/nelsayed/raw_data/ade/HPGL0242/processed/outputs/qsub/bt2idx-hpgl0242.qsubout -j oe -V -m n
echo "####Started /cbcb/lab/nelsayed/raw_data/ade/HPGL0242/processed/scripts/bt2idx-hpgl0242.sh at $(date)" >> outputs/log.txt
cd /cbcb/lab/nelsayed/raw_data/ade/HPGL0242/processed || exit

## Generating bowtie2 indexes for species: lamazonensis in /cbcbhomes/abelew/libraries/genome/indexes

if [ ! -r "/cbcbhomes/abelew/libraries/genome/lamazonensis.fa" ]; then
  ln -s /cbcbhomes/abelew/libraries/genome/lamazonensis.fasta /cbcbhomes/abelew/libraries/genome/lamazonensis.fa
fi
if [ ! -r "/cbcbhomes/abelew/libraries/genome/indexes/lamazonensis.fa" ]; then
  ln -s /cbcbhomes/abelew/libraries/genome/lamazonensis.fasta /cbcbhomes/abelew/libraries/genome/indexes/lamazonensis.fa
fi

bowtie2-build /cbcbhomes/abelew/libraries/genome/lamazonensis.fasta /cbcbhomes/abelew/libraries/genome/indexes/lamazonensis

## The following lines give status codes and some logging
echo $? > outputs/status/bt2idx-hpgl0242.status
echo "###Finished ${PBS_JODID} bt2idx-hpgl0242.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


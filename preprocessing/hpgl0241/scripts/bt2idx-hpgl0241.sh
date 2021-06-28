#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/old/hpgl0241
#PBS -N bt2idx-hpgl0241 -l mem=6gb -l walltime=10:00:00 -l ncpus=4
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/old/hpgl0241/outputs/qsub/bt2idx-hpgl0241.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/old/hpgl0241/scripts/bt2idx-hpgl0241.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/old/hpgl0241 || exit

## Generating bowtie2 indexes for species: ldonovani in /cbcbhomes/abelew/libraries/genome/indexes

if [ ! -r "/cbcbhomes/abelew/libraries/genome/ldonovani.fa" ]; then
  ln -s /cbcbhomes/abelew/libraries/genome/ldonovani.fasta /cbcbhomes/abelew/libraries/genome/ldonovani.fa
fi
if [ ! -r "/cbcbhomes/abelew/libraries/genome/indexes/ldonovani.fa" ]; then
  ln -s /cbcbhomes/abelew/libraries/genome/ldonovani.fasta /cbcbhomes/abelew/libraries/genome/indexes/ldonovani.fa
fi

bowtie2-build /cbcbhomes/abelew/libraries/genome/ldonovani.fasta /cbcbhomes/abelew/libraries/genome/indexes/ldonovani

## The following lines give status codes and some logging
echo $? > outputs/status/bt2idx-hpgl0241.status
echo "###Finished ${PBS_JODID} bt2idx-hpgl0241.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0658
#PBS -N hts_hsapiens-hpgl0658 -l mem=6gb -l walltime=10:00:00 -l ncpus=4
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0658/outputs/qsub/hts_hsapiens-hpgl0658.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0658/scripts/hts_hsapiens-hpgl0658.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0658 || exit

## Counting the number of hits in outputs/tophat_hsapiens/accepted_hits.bam for each feature found in /cbcbhomes/abelew/libraries/genome/hsapiens.gtf
## Is this stranded? no.  The defaults of htseq are:
##  --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union 

htseq-count -q -f bam -s no   \
  outputs/tophat_hsapiens/accepted_hits.bam /cbcbhomes/abelew/libraries/genome/hsapiens.gtf \
  1>outputs/tophat_hsapiens/accepted_hits.count 2>outputs/tophat_hsapiens/accepted_hits.error && \
    xz -9e outputs/tophat_hsapiens/accepted_hits.count

## The following lines give status codes and some logging
echo $? > outputs/status/hts_hsapiens-hpgl0658.status
echo "###Finished ${PBS_JODID} hts_hsapiens-hpgl0658.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


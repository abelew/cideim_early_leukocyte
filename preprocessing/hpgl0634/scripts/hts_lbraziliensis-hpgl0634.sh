#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0634
#PBS -N hts_lbraziliensis-hpgl0634 -l mem=6gb -l walltime=10:00:00 -l ncpus=4
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0634/outputs/qsub/hts_lbraziliensis-hpgl0634.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0634/scripts/hts_lbraziliensis-hpgl0634.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0634 || exit

## Counting the number of hits in outputs/tophat_lbraziliensis/accepted_hits.bam for each feature found in /cbcbhomes/abelew/libraries/genome/lbraziliensis.gff
## Is this stranded? no.  The defaults of htseq are:
##  --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union 

htseq-count -q -f bam -s no  -i ID  \
  outputs/tophat_lbraziliensis/accepted_hits.bam /cbcbhomes/abelew/libraries/genome/lbraziliensis.gff \
  1>outputs/tophat_lbraziliensis/accepted_hits.count 2>outputs/tophat_lbraziliensis/accepted_hits.error && \
    xz -9e outputs/tophat_lbraziliensis/accepted_hits.count

## The following lines give status codes and some logging
echo $? > outputs/status/hts_lbraziliensis-hpgl0634.status
echo "###Finished ${PBS_JODID} hts_lbraziliensis-hpgl0634.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


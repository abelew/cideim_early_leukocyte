#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0632
#PBS -N hts_lmexicana-hpgl0632 -l mem=6gb -l walltime=10:00:00 -l ncpus=4
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0632/outputs/qsub/hts_lmexicana-hpgl0632.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0632/scripts/hts_lmexicana-hpgl0632.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0632 || exit

## Counting the number of hits in outputs/tophat_lmexicana/accepted_hits.bam for each feature found in /cbcbhomes/abelew/libraries/genome/lmexicana.gff
## Is this stranded? no.  The attribute type (3rd column of the gff file) is exon
## and the comment field (10th column) used to name the feature is ID
htseq-count -q -f bam -s no -t exon -i ID \
  outputs/tophat_lmexicana/accepted_hits.bam /cbcbhomes/abelew/libraries/genome/lmexicana.gff \
  1>outputs/tophat_lmexicana/accepted_hits.count 2>outputs/tophat_lmexicana/accepted_hits.error && \
 pxz outputs/tophat_lmexicana/accepted_hits.count

## The following lines give status codes and some logging
echo $? > outputs/status/hts_lmexicana-hpgl0632.status
echo "###Finished ${PBS_JODID} hts_lmexicana-hpgl0632.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


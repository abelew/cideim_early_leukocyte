#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0633
#PBS -N th_lbraziliensis-hpgl0633 -l mem=8gb -l walltime=18:00:00 -l ncpus=4
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0633/outputs/qsub/th_lbraziliensis-hpgl0633.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0633/scripts/th_lbraziliensis-hpgl0633.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0633 || exit

## I still have no clue what I am doing when I use tophat...
## However, I know that -g 1 will allow only 1 hit in the case of multihits, but randomly place it
## From the manual:  "If there are more alignments with the same score than this
## number, TopHat will randomly report only this many alignments"
## -N 1 will discard anything with >1 mismatch (default is 2)
## -r adjusts the allowable mean distance between the paired reads
## --mate-std-dev sets the deviation of -r
## --microexon-search will tell it to search short exons for reads >=50

mkdir -p outputs/tophat_lbraziliensis && tophat  -g 1  \
  -G /cbcbhomes/abelew/libraries/genome/lbraziliensis.gff \
  --b2-very-sensitive -p 4 -o outputs/tophat_lbraziliensis \
/cbcbhomes/abelew/libraries/genome/indexes/lbraziliensis \
  hpgl0633_forward-trimmed.fastq.gz hpgl0633_reverse-trimmed.fastq.gz && \
  samtools sort -l 9 -n outputs/tophat_lbraziliensis/accepted_hits.bam outputs/tophat_lbraziliensis/accepted_sorted && \
  mv outputs/tophat_lbraziliensis/accepted_sorted.bam outputs/tophat_lbraziliensis/accepted_hits.bam && \
  samtools index outputs/tophat_lbraziliensis/accepted_hits.bam && \
  samtools sort -l 9 -n outputs/tophat_lbraziliensis/unmapped.bam outputs/tophat_lbraziliensis/unmapped_sorted && \
  mv outputs/tophat_lbraziliensis/unmapped_sorted.bam outputs/tophat_lbraziliensis/unmapped.bam && \
  samtools index outputs/tophat_lbraziliensis/unmapped.bam

## The following lines give status codes and some logging
echo $? > outputs/status/th_lbraziliensis-hpgl0633.status
echo "###Finished ${PBS_JODID} th_lbraziliensis-hpgl0633.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


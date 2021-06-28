#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/old/hpgl0244
#PBS -N trim-hpgl0244 -l mem=6gb -l walltime=4:00:00 -l ncpus=4
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/old/hpgl0244/outputs/qsub/trim-hpgl0244.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/old/hpgl0244/scripts/trim-hpgl0244.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/old/hpgl0244 || exit

## This call to trimomatic removes illumina and epicentre adapters from hpgl0244_forward.fastq.gz:hpgl0244_reverse.fastq.gz.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.

## Trimomatic_Pairwise: In case a trimming needs to be redone...
if [[ ! -r "hpgl0244_forward.fastq.gz" ]]; then
  if [[ -r "sequences/hpgl0244_forward.fastq.xz" ]]; then
    mv sequences/hpgl0244_forward.fastq.xz . && pxz -d hpgl0244_forward.fastq.xz && pigz hpgl0244_forward.fastq && mv sequences/hpgl0244_reverse.fastq.xz . && pxz -d hpgl0244_reverse.fastq.xz && pigz hpgl0244_reverse.fastq
  else
    echo "Missing files. Did not find hpgl0244_forward.fastq.gz nor sequences/hpgl0244_forward.fastq.xz"
    exit 1
  fi
fi
trimomatic PE -threads 1 -phred33 hpgl0244_forward.fastq.gz hpgl0244_reverse.fastq.gz hpgl0244_forward-trimmed_paired.fastq hpgl0244_forward-trimmed_unpaired.fastq hpgl0244_reverse-trimmed_paired.fastq hpgl0244_reverse-trimmed_unpaired.fastq ILLUMINACLIP:/cbcbhomes/abelew/libraries/adapters.fa:2:20:4 SLIDINGWINDOW:4:25 1>outputs/hpgl0244-trimomatic.out 2>&1
excepted=$(grep "Exception" outputs/hpgl0244-trimomatic.out)
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "${excepted}" != "" ]]; then
  trimomatic PE -threads 1 -phred33 hpgl0244_forward.fastq.gz hpgl0244_reverse.fastq.gz hpgl0244_forward-trimmed_paired.fastq hpgl0244_forward-trimmed_unpaired.fastq hpgl0244_reverse-trimmed_paired.fastq hpgl0244_reverse-trimmed_unpaired.fastq SLIDINGWINDOW:4:25 1>>outputs/hpgl0244-trimomatic.out 2>&1
fi
sleep 10
mv hpgl0244_forward-trimmed_paired.fastq hpgl0244_forward-trimmed.fastq && mv hpgl0244_reverse-trimmed_paired.fastq hpgl0244_reverse-trimmed.fastq

## The following lines give status codes and some logging
echo $? > outputs/status/trim-hpgl0244.status
echo "###Finished ${PBS_JODID} trim-hpgl0244.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt

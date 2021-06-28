#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0662/processed
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=trim_hpgl0662_forward
#SBATCH --mem=40G
#SBATCH --cpus-per-task=3
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0662/processed/outputs/trim_hpgl0662_forward.sbatchout


echo "## Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0662/processed/scripts/05trim_hpgl0662_forward.sh at $(date) on $(hostname) with id ${SLURM_JOBID}." >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0662/processed || exit

## This call to trimomatic removes illumina and epicentre adapters from hpgl0662_forward.fastq.xz:hpgl0662_reverse.fastq.xz.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.

## Trimomatic_Pairwise: In case a trimming needs to be redone...
if [[ ! -r "hpgl0662_forward.fastq.xz" ]]; then
  if [[ -r "sequences/hpgl0662_forward.fastq.xz" ]]; then
    mv sequences/hpgl0662_forward.fastq.xz . && pxz -d hpgl0662_forward.fastq.xz && pigz hpgl0662_forward.fastq &&\
       mv sequences/hpgl0662_reverse.fastq.xz . && pxz -d hpgl0662_reverse.fastq.xz && pigz hpgl0662_reverse.fastq
  else
    echo "Missing files. Did not find hpgl0662_forward.fastq.xz nor sequences/hpgl0662_forward.fastq.xz"
    exit 1
  fi
fi
## Note that trimomatic prints all output and errors to STDERR, so send both to output
trimmomatic PE \
  -threads 1 \
  -phred33 \
  <(less hpgl0662_forward.fastq.xz) <(less hpgl0662_reverse.fastq.xz) \
  hpgl0662_forward-trimmed_paired.fastq.gz hpgl0662_forward-trimmed_unpaired.fastq.gz \
  hpgl0662_reverse-trimmed_paired.fastq.gz hpgl0662_reverse-trimmed_unpaired.fastq.gz \
   ILLUMINACLIP:/fs/cbcb-software/RedHat-7-x86_64/common/local/perl/common/5.28.1/lib/perl5/x86_64-linux/auto/Bio/Adventure/genome/adapters.fa:2:30:10:2:keepBothReads \
  SLIDINGWINDOW:4:25 MINLEN:40 \
  1>outputs/hpgl0662-trimomatic.out 2>&1
excepted=$(grep "Exception" outputs/hpgl0662-trimomatic.out)
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "${excepted}" != "" ]]; then
  trimmomatic PE \
    -threads 1 \
    -phred33 \
    <(less hpgl0662_forward.fastq.xz) <(less hpgl0662_reverse.fastq.xz) \
    hpgl0662_forward-trimmed_paired.fastq.gz hpgl0662_forward-trimmed_unpaired.fastq.gz \
    hpgl0662_reverse-trimmed_paired.fastq.gz hpgl0662_reverse-trimmed_unpaired.fastq.gz \
     SLIDINGWINDOW:4:25 MINLEN:50\
    1>outputs/hpgl0662-trimomatic.out 2>&1
fi
sleep 10
mv hpgl0662_forward-trimmed_paired.fastq.gz hpgl0662_forward-trimmed.fastq.gz && mv hpgl0662_reverse-trimmed_paired.fastq.gz hpgl0662_reverse-trimmed.fastq.gz
ln -s hpgl0662_forward-trimmed.fastq.gz r1_trimmed.fastq.gz
ln -s hpgl0662_reverse-trimmed.fastq.gz r2_trimmed.fastq.gz

## The following lines give status codes and some logging
echo $? > outputs/status/trim_hpgl0662_forward.status
echo "## Finished ${SLURM_JOBID} 05trim_hpgl0662_forward.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt


walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


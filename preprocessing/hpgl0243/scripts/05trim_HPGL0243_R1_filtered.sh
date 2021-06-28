#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0243/processed
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=trim_HPGL0243_R1_filtered
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0243/processed/outputs/trim_HPGL0243_R1_filtered.sbatchout

echo "####Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0243/processed/scripts/05trim_HPGL0243_R1_filtered.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0243/processed || exit

## This call to trimomatic removes illumina and epicentre adapters from HPGL0243_R1_filtered.fastq.gz:HPGL0243_R2_filtered.fastq.gz.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.

## Trimomatic_Pairwise: In case a trimming needs to be redone...
if [[ ! -r "HPGL0243_R1_filtered.fastq.gz" ]]; then
  if [[ -r "sequences/HPGL0243_R1_filtered.fastq.xz" ]]; then
    mv sequences/HPGL0243_R1_filtered.fastq.xz . && pxz -d HPGL0243_R1_filtered.fastq.xz && pigz HPGL0243_R1_filtered.fastq &&\
       mv sequences/HPGL0243_R2_filtered.fastq.xz . && pxz -d HPGL0243_R2_filtered.fastq.xz && pigz HPGL0243_R2_filtered.fastq
  else
    echo "Missing files. Did not find HPGL0243_R1_filtered.fastq.gz nor sequences/HPGL0243_R1_filtered.fastq.xz"
    exit 1
  fi
fi
## Note that trimomatic prints all output and errors to STDERR, so send both to output
trimmomatic PE -threads 1 -phred33 HPGL0243_R1_filtered.fastq.gz HPGL0243_R2_filtered.fastq.gz HPGL0243_R1_filtered-trimmed_paired.fastq.gz HPGL0243_R1_filtered-trimmed_unpaired.fastq.gz HPGL0243_R2_filtered-trimmed_paired.fastq.gz HPGL0243_R2_filtered-trimmed_unpaired.fastq.gz \
    ILLUMINACLIP:/cbcb/sw/RedHat-7-x86_64/common/local/perl/common/5.26/lib/auto/share/module/Bio-Adventure/genome/adapters.fa:2:20:4 SLIDINGWINDOW:4:25 \
    1>outputs/HPGL0243_R1_filtered-trimomatic.out 2>&1
excepted=$(grep "Exception" outputs/HPGL0243_R1_filtered-trimomatic.out)
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "${excepted}" != "" ]]; then
  trimmomatic PE -threads 1 -phred33 HPGL0243_R1_filtered.fastq.gz HPGL0243_R2_filtered.fastq.gz HPGL0243_R1_filtered-trimmed_paired.fastq.gz HPGL0243_R1_filtered-trimmed_unpaired.fastq.gz HPGL0243_R2_filtered-trimmed_paired.fastq.gz HPGL0243_R2_filtered-trimmed_unpaired.fastq.gz SLIDINGWINDOW:4:25 \
    1>outputs/HPGL0243_R1_filtered-trimomatic.out 2>&1
fi
sleep 10
mv HPGL0243_R1_filtered-trimmed_paired.fastq.gz HPGL0243_R1_filtered-trimmed.fastq.gz && mv HPGL0243_R2_filtered-trimmed_paired.fastq.gz HPGL0243_R2_filtered-trimmed.fastq.gz

## The following lines give status codes and some logging
echo $? > outputs/status/trim_HPGL0243_R1_filtered.status
echo "###Finished ${SLURM_JOBID} 05trim_HPGL0243_R1_filtered.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0662/processed
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=sal_hg38_91
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0662/processed/outputs/sal_hg38_91.sbatchout


echo "####Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0662/processed/scripts/30sal_hg38_91.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0662/processed || exit

## This is a salmon pseudoalignment of forward.fastq.gz:reverse.fastq.gz against
## /cbcbhomes/abelew/libraries/genome/indexes/hg38_91_salmon_index.
## This jobs depended on: 

mkdir -p outputs/salmon_hg38_91 && sleep 3 && \
salmon quant -i /cbcbhomes/abelew/libraries/genome/indexes/hg38_91_salmon_index \
  -l A --gcBias  \
   -1 <(less forward.fastq.gz) -2 <(less reverse.fastq.gz)  \
  -o outputs/salmon_hg38_91 \
  2>outputs/salmon_hg38_91/salmon.err 1>outputs/salmon_hg38_91/salmon.out

## The following lines give status codes and some logging
echo $? > outputs/status/sal_hg38_91.status
echo "###Finished ${SLURM_JOBID} 30sal_hg38_91.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


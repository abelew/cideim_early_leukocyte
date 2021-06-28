#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0242/processed
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=hisat2_linfantum_v46_r1
#SBATCH --mem=45G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0242/processed/outputs/hisat2_linfantum_v46_r1.sbatchout


echo "####Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0242/processed/scripts/15hisat2_linfantum_v46_r1.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0242/processed || exit

## This is a hisat2 alignment of  -1 <(less /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0242/processed/r1.fastq.xz) -2 <(less /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0242/processed/r2.fastq.xz)  against
## /cbcbhomes/abelew/libraries/genome/indexes/linfantum_v46 using arguments: .
## This jobs depended on: 

mkdir -p outputs/hisat2_linfantum_v46 && \
  sleep 3 && \
  hisat2 -x /cbcbhomes/abelew/libraries/genome/indexes/linfantum_v46  \
    -p 4 \
    -q   -1 <(less /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0242/processed/r1.fastq.xz) -2 <(less /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0242/processed/r2.fastq.xz)  \
    --phred33 \
    --un outputs/hisat2_linfantum_v46/r1_unaligned_linfantum_v46.fastq \
    --al outputs/hisat2_linfantum_v46/r1_aligned_linfantum_v46.fastq \
    -S outputs/hisat2_linfantum_v46/r1.sam \
    2>outputs/hisat2_linfantum_v46/r1.err \
    1>outputs/hisat2_linfantum_v46/r1.out

## The following lines give status codes and some logging
echo $? > outputs/status/hisat2_linfantum_v46_r1.status
echo "###Finished ${SLURM_JOBID} 15hisat2_linfantum_v46_r1.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


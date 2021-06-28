#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0658/processed
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=bt2_lpanamensis_v36_r1
#SBATCH --mem=28G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0658/processed/outputs/bt2_lpanamensis_v36_r1.sbatchout


echo "## Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0658/processed/scripts/15bt2_lpanamensis_v36_r1.sh at $(date) on $(hostname) with id ${SLURM_JOBID}." >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0658/processed || exit

## This is a bowtie2 alignment of  -1 <(less /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0658/processed/r1.fastq.xz) -2 <(less /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0658/processed/r2.fastq.xz)  against
## /cbcbhomes/abelew/libraries/genome/indexes/lpanamensis_v36 using arguments:  --very-sensitive -L 14 .
## This jobs depended on: 

mkdir -p outputs/bowtie2_lpanamensis_v36 && \
  sleep 3 && \
  bowtie2 -x /cbcbhomes/abelew/libraries/genome/indexes/lpanamensis_v36  --very-sensitive -L 14  \
    -p 4 \
    -q   -1 <(less /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0658/processed/r1.fastq.xz) -2 <(less /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0658/processed/r2.fastq.xz)  \
    --un outputs/bowtie2_lpanamensis_v36/r1_unaligned_lpanamensis_v36.fastq \
    --al outputs/bowtie2_lpanamensis_v36/r1_aligned_lpanamensis_v36.fastq \
    -S outputs/bowtie2_lpanamensis_v36/r1.sam \
    2>outputs/bowtie2_lpanamensis_v36/r1.err \
    1>outputs/bowtie2_lpanamensis_v36/r1.out

## The following lines give status codes and some logging
echo $? > outputs/status/bt2_lpanamensis_v36_r1.status
echo "## Finished ${SLURM_JOBID} 15bt2_lpanamensis_v36_r1.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt


walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


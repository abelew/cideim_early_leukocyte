#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed
#SBATCH --partition=dpart
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --time=18:00:00
#SBATCH --job-name=xzun_bt2_r1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed/outputs/xzun_bt2_r1.sbatchout


echo "## Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed/scripts/16xzun_bt2_r1.sh at $(date) on $(hostname) with id ${SLURM_JOBID}." >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed || exit

## Compressing the sequences which failed to align against /cbcbhomes/abelew/libraries/genome/indexes/lpanamensis_v36 using options  --very-sensitive -L 14 

nice -n 20 xz -f -9e outputs/bowtie2_lpanamensis_v36/r1_unaligned_lpanamensis_v36.fastq   2>outputs/bowtie2_lpanamensis_v36/xzun_bt2_r1.err \
   1>outputs/bowtie2_lpanamensis_v36/xzun_bt2_r1.out
## The following lines give status codes and some logging
echo $? > outputs/status/xzun_bt2_r1.status
echo "## Finished ${SLURM_JOBID} 16xzun_bt2_r1.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt


walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


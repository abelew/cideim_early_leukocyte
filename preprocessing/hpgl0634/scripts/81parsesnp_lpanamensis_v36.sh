#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=parsesnp_lpanamensis_v36
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed/outputs/parsesnp_lpanamensis_v36.sbatchout


echo "## Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed/scripts/81parsesnp_lpanamensis_v36.sh at $(date) on $(hostname) with id ${SLURM_JOBID}." >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed || exit


## Parse the SNP data and generate a modified lpanamensis_v36 genome.
##  This should read the file:
## outputs/vcfutils_lpanamensis_v36/r1.bcf
##  and provide 4 new files:
## outputs/vcfutils_lpanamensis_v36/r1_lpanamensis_v36_count.txt
## outputs/vcfutils_lpanamensis_v36/r1_lpanamensis_v36_all.txt
## outputs/vcfutils_lpanamensis_v36/r1_lpanamensis_v36_pct.txt
## and a modified genome: outputs/vcfutils_lpanamensis_v36/r1_lpanamensis_v36_modified.fasta

/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0634/processed/scripts/81parsesnp_lpanamensis_v36.pl

## The following lines give status codes and some logging
echo $? > outputs/status/parsesnp_lpanamensis_v36.status
echo "## Finished ${SLURM_JOBID} 81parsesnp_lpanamensis_v36.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt


walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


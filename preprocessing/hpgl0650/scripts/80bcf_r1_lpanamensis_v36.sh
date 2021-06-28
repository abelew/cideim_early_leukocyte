#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0650/processed
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=bcf_r1_lpanamensis_v36
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0650/processed/outputs/bcf_r1_lpanamensis_v36.sbatchout


echo "## Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0650/processed/scripts/80bcf_r1_lpanamensis_v36.sh at $(date) on $(hostname) with id ${SLURM_JOBID}." >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0650/processed || exit

## Use samtools, bcftools, and vcfutils to get some idea about how many variant positions are in the data.
mkdir -p outputs/vcfutils_lpanamensis_v36
echo "Started samtools sort at $(date)" >> outputs/vcfutils_lpanamensis_v36.out
samtools sort -l 9 -@ 4 outputs/bowtie2_lpanamensis_v36/r1.bam -o outputs/vcfutils_lpanamensis_v36/r1.bam 2>outputs/samtools_sort.out 1>&2
if [ "$?" -ne "0" ]; then
    echo "samtools sort failed."
exit 1
fi

if [ ! -r "/cbcbhomes/abelew/libraries/genome/lpanamensis_v36.fasta.fai" ]; then
    samtools faidx /cbcbhomes/abelew/libraries/genome/lpanamensis_v36.fasta
fi
samtools mpileup -uvf /cbcbhomes/abelew/libraries/genome/lpanamensis_v36.fasta 2>samtools_mpileup.err \
    outputs/vcfutils_lpanamensis_v36/r1.bam |\
  bcftools call -c - 2>bcftools_call.err |\
  bcftools view -l 9 -o outputs/vcfutils_lpanamensis_v36/r1.bcf -O b - \
    2>outputs/vcfutils_lpanamensis_v36/r1_summary.err
if [ "$?" -ne "0" ]; then
    echo "mpileup/bcftools failed."
    exit 1
fi
bcftools index outputs/vcfutils_lpanamensis_v36/r1.bcf 2>bcftools_index.err
echo "Successfully finished." >> outputs/vcfutils_lpanamensis_v36.out

## The following lines give status codes and some logging
echo $? > outputs/status/bcf_r1_lpanamensis_v36.status
echo "## Finished ${SLURM_JOBID} 80bcf_r1_lpanamensis_v36.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt


walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


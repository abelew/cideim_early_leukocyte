#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0650/processed
#SBATCH --partition=dpart
#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=ht2st_hisat2_r1
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0650/processed/outputs/ht2st_hisat2_r1.sbatchout


echo "#### Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0650/processed/scripts/18ht2st_hisat2_r1.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0650/processed || exit

## This is a stupidly simple job to collect alignment statistics.

if [ ! -e "outputs/hisat2_stats.csv" ]; then
    echo "original reads, single hits, failed reads, multi-hits, rpm" > outputs/hisat2_stats.csv
fi
original_reads_tmp=$(grep " reads; of these" "outputs/hisat2_leishmania_rna_virus_1/r1.err" 2>/dev/null | awk '{print $1}' | sed 's/ //g')
original_reads=${original_reads_tmp:-0}
one_align_tmp=$(grep " aligned exactly 1 time" "outputs/hisat2_leishmania_rna_virus_1/r1.err" | awk '{print $1}' | sed 's/ .*//g')
one_align=${one_align_tmp:-0}
failed_tmp=$(grep " aligned 0 times" "outputs/hisat2_leishmania_rna_virus_1/r1.err" | awk '{print $1}' | sed 's/ .*//g')
failed=${failed_tmp:-0}
sampled_tmp=$(grep " aligned >1 times" "outputs/hisat2_leishmania_rna_virus_1/r1.err" | awk '{print $1}' | sed 's/ .*//g')
sampled=${sampled_tmp:-0}
rpm_tmp=$(perl -e "printf(1000000 / $(( ${one_align} + ${sampled} )) ) " 2>/dev/null)
rpm=${rpm_tmp:-0}
stat_string=$(printf "r1,%s,%s,%s,%s,%s" "${original_reads}" "${one_align}" "${failed}" "${sampled}" "${rpm}")
echo "$stat_string" >> outputs/hisat2_stats.csv
## The following lines give status codes and some logging
echo $? > outputs/status/ht2st_hisat2_r1.status
echo "### Finished ${SLURM_JOBID} 18ht2st_hisat2_r1.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


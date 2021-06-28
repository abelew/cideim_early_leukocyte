#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0655/processed
#SBATCH --partition=dpart
#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=sastats_hg38_99
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0655/processed/outputs/sastats_hg38_99.sbatchout


echo "####Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0655/processed/scripts/33sastats_hg38_99.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0655/processed || exit

## This is a stupidly simple job to collect salmon alignment statistics.

if [ ! -r "outputs/salmon_stats.csv" ]; then
  echo "basename,species,fragments,assigned,consistent,inconsistent,bias" > outputs/salmon_stats.csv
fi
reads_tmp=$(grep "^num_compatible" outputs/salmon_hg38_99/lib_format_counts.json | awk '{print $3}' | sed 's/^ *//g')
reads=${reads_tmp:-0}
aligned_tmp=$(grep "^num_assigned" outputs/salmon_hg38_99/lib_format_counts.json | awk '{print $3}' | sed 's/^ *//g')
aligned=${aligned_tmp:-0}
consistent_tmp=$(grep "^concordant" outputs/salmon_hg38_99/lib_format_counts.json | awk '{print $3}' | sed 's/^ *//g')
consistent=${consistent_tmp:-0}
inconsistent_tmp=$(grep "^inconsistent" outputs/salmon_hg38_99/lib_format_counts.json | awk '{print $3}' | sed 's/^ *//g')
inconsistent=${inconsistent_tmp:-0}
bias_tmp=$(grep "^mapping_bias" outputs/salmon_hg38_99/lib_format_counts.json | awk '{print $3}' | sed 's/^ *//g')
bias=${bias_tmp:-0}
stat_string=$(printf "r1,hg38_99,%s,%s,%s,%s,%s" "${reads}" "${aligned}" "${consistent}" "${inconsistent}" "${bias}")
echo "$stat_string" >> "outputs/salmon_stats.csv"
## The following lines give status codes and some logging
echo $? > outputs/status/sastats_hg38_99.status
echo "###Finished ${SLURM_JOBID} 33sastats_hg38_99.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


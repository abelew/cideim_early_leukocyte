#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0243/processed
#SBATCH --partition=dpart
#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=trst_HPGL0243_R1_filtered
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0243/processed/outputs/trst_HPGL0243_R1_filtered.sbatchout

echo "####Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0243/processed/scripts/06trst_HPGL0243_R1_filtered.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0243/processed || exit

## This is a stupidly simple job to collect trimomatic statistics

if [ ! -r outputs/trimomatic_stats.csv ]; then
  echo "total_reads,surviving_both,surviving_forward,surviving_reverse,dropped_reads" > outputs/trimomatic_stats.csv
fi
total_reads_tmp=$(grep "^Input Read Pairs" outputs/HPGL0243_R1_filtered-trimomatic.out | awk '{print $4}')
total_reads=${total_reads_tmp:-0}
surviving_both_tmp=$(grep "^Input Read Pairs" outputs/HPGL0243_R1_filtered-trimomatic.out | awk '{print $7}')
surviving_both=${surviving_both_tmp:-0}
surviving_forward_tmp=$(grep "^Input Read Pairs" outputs/HPGL0243_R1_filtered-trimomatic.out | awk '{print $12}')
surviving_forward=${surviving_forward_tmp:-0}
surviving_reverse_tmp=$(grep "^Input Read Pairs" outputs/HPGL0243_R1_filtered-trimomatic.out | awk '{print $17}')
surviving_reverse=${surviving_reverse_tmp:-0}
dropped_reads_tmp=$(grep "^Input Read Pairs" outputs/HPGL0243_R1_filtered-trimomatic.out | awk '{print $20}')
dropped_reads=${dropped_reads_tmp:-0}

stat_string=$(printf "HPGL0243_R1_filtered,%s,%s,%s,%s,%s" "${total_reads}" "${surviving_both}" "${surviving_forward}" "${surviving_reverse}" "${dropped_reads}")
echo "$stat_string" >> outputs/trimomatic_stats.csv

## The following lines give status codes and some logging
echo $? > outputs/status/trst_HPGL0243_R1_filtered.status
echo "###Finished ${SLURM_JOBID} 06trst_HPGL0243_R1_filtered.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


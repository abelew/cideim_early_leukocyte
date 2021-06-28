#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0245/processed
#SBATCH --partition=dpart
#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=htall_hisat2_r1
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0245/processed/outputs/htall_hisat2_r1.sbatchout


echo "####Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0245/processed/scripts/21htall_hisat2_r1.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0245/processed || exit

## Counting the number of hits in outputs/hisat2_ltropica590_v46/r1.bam for each feature found in /cbcbhomes/abelew/libraries/genome/ltropica590_v46.gff
## Is this stranded? no.  The defaults of htseq are:
##  --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union 


htseq-count  --help 2>&1 | tail -n 3
htseq-count \
  -q -f bam -s no  --idattr ID  --type exon \
  outputs/hisat2_ltropica590_v46/r1.bam \
  /cbcbhomes/abelew/libraries/genome/ltropica590_v46.gff \
  2>outputs/hisat2_ltropica590_v46/r1_htseq.err \
  1>outputs/hisat2_ltropica590_v46/r1.count && \
    xz -f -9e outputs/hisat2_ltropica590_v46/r1.count 2>outputs/hisat2_ltropica590_v46/r1_htseq.err.xz 1>outputs/hisat2_ltropica590_v46/r1.count.xz

## The following lines give status codes and some logging
echo $? > outputs/status/htall_hisat2_r1.status
echo "###Finished ${SLURM_JOBID} 21htall_hisat2_r1.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


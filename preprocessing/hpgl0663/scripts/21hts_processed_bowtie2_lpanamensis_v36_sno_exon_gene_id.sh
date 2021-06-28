#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0663/processed
#SBATCH --partition=dpart
#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=hts_processed_bowtie2_lpanamensis_v36_sno_exon_gene_id
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0663/processed/outputs/hts_processed_bowtie2_lpanamensis_v36_sno_exon_gene_id.sbatchout


echo "## Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0663/processed/scripts/21hts_processed_bowtie2_lpanamensis_v36_sno_exon_gene_id.sh at $(date) on $(hostname) with id ${SLURM_JOBID}." >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0663/processed || exit

## Counting the number of hits in outputs/bowtie2_lpanamensis_v36/r1.bam for each feature found in /cbcbhomes/abelew/libraries/genome/lpanamensis_v36.gff
## Is this stranded? no.  The defaults of htseq are:
##  --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union 


htseq-count  --help 2>&1 | tail -n 3
htseq-count \
  -q -f bam -s no  --type exon  --idattr gene_id \
  outputs/bowtie2_lpanamensis_v36/r1.bam \
  /cbcbhomes/abelew/libraries/genome/lpanamensis_v36.gff \
  2>outputs/bowtie2_lpanamensis_v36/r1_htseq.err \
  1>outputs/bowtie2_lpanamensis_v36/r1.count_lpanamensis_v36_sno_exon_gene_id.count && \
    xz -f -9e outputs/bowtie2_lpanamensis_v36/r1.count_lpanamensis_v36_sno_exon_gene_id.count 2>outputs/bowtie2_lpanamensis_v36/r1_htseq.err.xz 1>outputs/bowtie2_lpanamensis_v36/r1.count_lpanamensis_v36_sno_exon_gene_id.count.xz

## The following lines give status codes and some logging
echo $? > outputs/status/hts_processed_bowtie2_lpanamensis_v36_sno_exon_gene_id.status
echo "## Finished ${SLURM_JOBID} 21hts_processed_bowtie2_lpanamensis_v36_sno_exon_gene_id.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt


walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


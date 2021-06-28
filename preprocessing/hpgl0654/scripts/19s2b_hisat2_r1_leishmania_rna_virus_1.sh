#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0654/processed
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=s2b_hisat2_r1_leishmania_rna_virus_1
#SBATCH --mem=18G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0654/processed/outputs/s2b_hisat2_r1_leishmania_rna_virus_1.sbatchout


echo "#### Started /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0654/processed/scripts/19s2b_hisat2_r1_leishmania_rna_virus_1.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /fs/cbcb-lab/nelsayed/raw_data/ade/HPGL0654/processed || exit

## Converting the text sam to a compressed, sorted, indexed bamfile.
## Also printing alignment statistics to outputs/hisat2_leishmania_rna_virus_1/r1.bam.stats
## This job depended on: 1775720

if $(test ! -r outputs/hisat2_leishmania_rna_virus_1/r1.sam); then
    echo "Could not find the samtools input file."
    exit 1
fi
samtools view -u -t /cbcbhomes/abelew/libraries/genome/leishmania_rna_virus_1.fasta \
  -S outputs/hisat2_leishmania_rna_virus_1/r1.sam -o outputs/hisat2_leishmania_rna_virus_1/r1.bam \
  2>outputs/hisat2_leishmania_rna_virus_1/r1.bam.err 1>outputs/hisat2_leishmania_rna_virus_1/r1.bam.out && \
  samtools sort -l 9 outputs/hisat2_leishmania_rna_virus_1/r1.bam -o outputs/hisat2_leishmania_rna_virus_1/r1-sorted.bam \
  2>outputs/hisat2_leishmania_rna_virus_1/r1-sorted.err 1>outputs/hisat2_leishmania_rna_virus_1/r1-sorted.out && \
  rm outputs/hisat2_leishmania_rna_virus_1/r1.bam && \
  rm outputs/hisat2_leishmania_rna_virus_1/r1.sam && \
  mv outputs/hisat2_leishmania_rna_virus_1/r1-sorted.bam outputs/hisat2_leishmania_rna_virus_1/r1.bam && \
  samtools index outputs/hisat2_leishmania_rna_virus_1/r1.bam
## The following will fail if this is single-ended.
samtools view -b -f 2 -o outputs/hisat2_leishmania_rna_virus_1/r1-paired.bam outputs/hisat2_leishmania_rna_virus_1/r1.bam && \
  samtools index outputs/hisat2_leishmania_rna_virus_1/r1-paired.bam
bamtools stats -in outputs/hisat2_leishmania_rna_virus_1/r1.bam 2>outputs/hisat2_leishmania_rna_virus_1/r1.bam.stats 1>&2 && \
  bamtools stats -in outputs/hisat2_leishmania_rna_virus_1/r1-paired.bam 2>outputs/hisat2_leishmania_rna_virus_1/r1-paired.stats 1>&2
##bamtools filter -tag XM:0 -in outputs/hisat2_leishmania_rna_virus_1/r1.bam -out outputs/hisat2_leishmania_rna_virus_1/r1-sorted_nomismatch.bam &&
##  samtools index outputs/hisat2_leishmania_rna_virus_1/r1-sorted_nomismatch.bam

## The following lines give status codes and some logging
echo $? > outputs/status/s2b_hisat2_r1_leishmania_rna_virus_1.status
echo "### Finished ${SLURM_JOBID} 19s2b_hisat2_r1_leishmania_rna_virus_1.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt


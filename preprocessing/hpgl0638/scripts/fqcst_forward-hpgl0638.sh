#PBS -V -S /usr/bin/bash -q throughput
#PBS -d /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0638
#PBS -N fqcst_forward-hpgl0638 -l mem=1gb -l walltime=00:10:00 -l ncpus=1
#PBS -o /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0638/outputs/qsub/fqcst_forward-hpgl0638.qsubout -j oe -V -m n
echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0638/scripts/fqcst_forward-hpgl0638.sh at $(date)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq_lpanamensis/preprocessing/hpgl0638 || exit

## This is a stupidly simple job to collect alignment statistics.

if [ ! -r outputs/fastqc_forward_stats.csv ]; then
  echo "total_reads,poor_quality,per_quality,per_base_content,per_sequence_gc,per_base_n,per_seq_length,over_rep,adapter_content,kmer_content" > outputs/fastqc_forward_stats.csv
fi
total_reads_tmp=$(grep "^Total Sequences" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
total_reads=${total_reads_tmp:-0}
poor_quality_tmp=$(grep "^Sequences flagged as poor quality" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
poor_quality=${poor_quality_tmp:-0}
per_quality_tmp=$(grep "Per base sequence quality" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
per_quality=${per_quality_tmp:-0}
per_base_content_tmp=$(grep "Per base sequence content" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
per_base_content=${per_base_content_tmp:-0}
per_sequence_gc_tmp=$(grep "Per sequence GC content" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
per_sequence_gc=${per_sequence_gc_tmp:-0}
per_base_n_tmp=$(grep "Per base N content" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
per_base_n=${per_base_n_tmp:-0}
per_seq_length_tmp=$(grep "Sequence Length Distribution" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
per_seq_length=${per_seq_length_tmp:-0}
over_rep_tmp=$(grep "Overrepresented sequences" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
over_rep=${over_rep_tmp:-0}
adapter_content_tmp=$(grep "Adapter Content" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
adapter_content=${adapter_content_tmp:-0}
kmer_content_tmp=$(grep "Kmer Content" outputs/fastqc/hpgl0638_forward.gz_fastqc/fastqc_data.txt | awk -F '\\t' '{print $2}')
kmer_content=${kmer_content_tmp:-0}

stat_string=$(printf "hpgl0638,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" "${total_reads}" "${poor_quality}" "${per_quality}" "${per_base_content}" "${per_sequence_gc}" "${per_base_n}" "${per_seq_length}" "${over_rep}" "${adapter_content}" "${kmer_content}")
echo "$stat_string" >> outputs/fastqc_forward_stats.csv

## The following lines give status codes and some logging
echo $? > outputs/status/fqcst_forward-hpgl0638.status
echo "###Finished ${PBS_JODID} fqcst_forward-hpgl0638.sh at $(date), it took $(( $SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
echo "####PBS walltime used by ${PBS_JOBID} was: ${walltime:-null}" >> outputs/log.txt
mem=$(qstat -f -t | grep "${PBS_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
echo "####PBS memory used by ${PBS_JOBID} was: ${mem:-null}" >> outputs/log.txt
vmmemory=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
echo "####PBS vmemory used by ${PBS_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
cputime=$(qstat -f -t "${PBS_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
echo "####PBS cputime used by ${PBS_JOBID} was: ${cputime:-null}" >> outputs/log.txt
##qstat -f -t ${PBS_JOBID} >> outputs/log.txt


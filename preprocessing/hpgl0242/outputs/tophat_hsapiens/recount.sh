#!/usr/bin/env bash

htseq-count -q -f bam -s no accepted_paired.bam ~/scratch/libraries/genome/hsapiens.gtf > accepted_paired.count && xz -f -9e accepted_paired.count


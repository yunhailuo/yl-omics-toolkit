#!/bin/bash

# Get the full directory name of the script no matter where it is being called from.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-which-directory-it-is-stored-in/246128#246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python $DIR/../../yl_fastq_tools.py \
       -f input.fastq \
       -s 18: \
       -r 1:5,12:17 \
       -i 6:11 \
       -b CGTGAT,ACATCG,CACTGT:Sample1_processed.fastq\;CTGATC:Sample2_processed.fastq\;AAGCTA,GCGGAC:Sample3_processed.fastq

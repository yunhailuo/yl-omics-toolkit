#!/bin/bash

# Get the full directory name of the script no matter where it is being called from.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-which-directory-it-is-stored-in/246128#246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change PICARD to your picard.jar location
PICARD=path/picard.jar

python $DIR/../../yl_fastq_tools.py \
       -f $DIR/input.fastq \
       -s 31:50 \
       -r 1:9 \
       -o $DIR/input_randomer.fastq

bowtie2-build reference.fasta reference
bowtie2 -L 20 -N 0 --no-1mm-upfront \
        -x reference \
        -U $DIR/input_randomer.fastq \
        -S $DIR/input.sam
rm *.bt2

python $DIR/../../yl_sam_tools.py \
       -s $DIR/input.sam \
       -o $DIR/input_bcqt.sam

java -Xmx2g -jar $PICARD SortSam I=$DIR/input_bcqt.sam \
                                 O=$DIR/input_bcqt_sorted.sam \
                                 SORT_ORDER=queryname

java -Xmx2g -jar $PICARD MarkDuplicates I=$DIR/input_bcqt_sorted.sam \
                                        O=$DIR/input_bcqt_sorted_dedup.sam \
                                        M=$DIR/picard_dedup.stats \
                                        BARCODE_TAG=BC \
                                        REMOVE_DUPLICATES=True

#!/bin/bash

# Get the full directory name of the script no matter where it is being called from.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-which-directory-it-is-stored-in/246128#246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change PICARD to your picard.jar location
PICARD=path/picard.jar

# Extract randomer from fastq
python $DIR/../../yl_fastq_tools.py \
       -f $DIR/input.fastq \
       -s 31:50 \
       -r 1:9 \
       -o $DIR/input_randomer.fastq
# Build bowtie2 index for alignment
bowtie2-build reference.fasta reference
# Bowtie2 alignment
bowtie2 -L 20 -N 0 --no-1mm-upfront \
        -x reference \
        -U $DIR/input_randomer.fastq \
        -S $DIR/input.sam
# For this example, remove bowtie2 index to save space
rm *.bt2

# Attach randomer to the BC&QT tag in the sam file
python $DIR/../../yl_sam_tools.py \
       -s $DIR/input.sam \
       -o $DIR/input_bcqt.sam

# Picard MarkDuplicates takes "either coordinate-sorted or query-sorted inputs."
# "However, when the input is query-sorted (actually query-grouped), then "
# "unmapped mates and secondary/supplementary reads are not excluded from the "
# "duplication test and can be marked as duplicate reads." It doesn't matter for
# this example though.
java -Xmx2g -jar $PICARD SortSam I=$DIR/input_bcqt.sam \
                                 O=$DIR/input_bcqt_sorted.sam \
                                 SORT_ORDER=queryname
# Picard MarkDuplicates with BARCODE_TAG option
java -Xmx2g -jar $PICARD MarkDuplicates I=$DIR/input_bcqt_sorted.sam \
                                        O=$DIR/input_bcqt_sorted_dedup.sam \
                                        M=$DIR/picard_dedup.stats \
                                        BARCODE_TAG=BC \
                                        REMOVE_DUPLICATES=True

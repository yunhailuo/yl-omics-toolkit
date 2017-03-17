# yl-omics-toolkit
A collection of bioinformatic tools customized to my needs.

**Table of Contents**
- [Usage](#user-content-usage)
  - [yl_fastq_tools](#user-content-yl_fastq_tools)
  - [yl_sam_tools](#user-content-yl_sam_tools)
- [Process randomer for Picard on marking PCR duplicates](#user-content-process-randomer-for-picard-on-marking-pcr-duplicates)
  - [1. Extract randomer from fastq and align to reference](#user-content-1-extract-randomer-from-fastq-and-align-to-reference)
  - [2. Attach randomer to the BC&QT tag in the sam file](#user-content-2-attach-randomer-to-the-bcqt-tag-in-the-sam-file)
  - [3. Picard MarkDuplicates with BARCODE_TAG option](#user-content-3-picard-markduplicates-with-barcode_tag-option)

## Usage
### yl_fastq_tools
`yl_fastq_tools.py [-h] -f FASTQ [-s SEQUENCE] [-i BCINDEX] [-b BCFILE] [-r RANDOMER] [-o OUTPUT]`
### Arguments
&nbsp;&nbsp;&nbsp;Argument&nbsp;&nbsp;&nbsp;|Help message
---|---
-f, --fastq|One fastq file from Illumina sequencing. Fastq files coming from other platforms have not been tested. Though pair end reads can be processed separately, it is not recommended since it may uncouple the pair.
-s, --sequence|Basepair position for target sequences which should be kept for downstream analysis. Support both continuous and gapped regions. Use syntax 'start1:stop1,start2:stop2,...'. Both the start and the end are included. If a stop is omitted, start-stop pairs after it will be ignored and sequencing bases will be kept up to the end. If this option is omitted, the whole sequence will be kept. However, at least one of the '-s', '-i'+'-b' and '-r' options must be provided.
-i, --bcindex|Basepair position for library barcodes which help demultiplex reads. Barcode sequences and corresponding output file names must be provided in the '-b' option (check it for details). Support both continuous and gapped regions. Use syntax 'start1:stop1,start2:stop2,...'. The syntax specification is the same as '-s' option above (check it for details). At least one of the '-s', '-i'+'-b' and '-r' options must be provided.
-b, --bcfile|Barcode sequences and corresponding output file names for demultiplexing. Use syntax 'barcode1-1,barcode1-2,barcode1-3,...:filename1;barcode2-1,barcode2-2,barcode2-3,...:filename2;...'. No space is allowed. At least one of the '-s', '-i'+'-b' and '-r' options must be provided.
-r, --randomer|Basepair position for randomer (molecular barcode) which help remove PCR duplicates generated during library preparation. It will be extracted and saved in corresponding reads' name temporarily. Later in the pipeline (in the sam file), it can be moved to the BC&QT tag which can be utilized by the Picard tool for marking/removing duplicates. Support both continuous and gapped regions. Use syntax 'start1:stop1,start2:stop2,...'. The syntax specifications are the same as '-s' option above (check it for details). At least one of the '-s', '-i'+'-b' and '-r' options must be provided.
-o, --output|Default output file. If '-i' and '-b' options are provided, the file specified by this option will be used for collecting reads failed to matchany of the barcodes. When this option is omitted, unmatched reads will be outputted to a fastq file named with input fastq file name plus a "_index_undetermined.fastq" postfix. If '-i' and '-b' options are omitted, all processed reads will be outputted to one single fastq file specified by this option. When this optioned is omitted, all processed reads will be outputted to a fastq file named with input fastq file name plus a "_processed.fastq" postfix.
-h, --help|Show this help message and exit.
### yl_sam_tools
`yl_sam_tools.py [-h] -s SAM [-o OUTSAM]`
### Arguments
&nbsp;&nbsp;&nbsp;Argument&nbsp;&nbsp;&nbsp;|Help message
---|---
-s --sam|One aligned sam file which needs to attach barcodes at BC&QT tags. The barcode and its sequencing quality should have been stored in the read's name, separated by ":". Check "yl_fastq_tools" module for detailed formats.
-o, --outsam|Output unsorted sam file. Default is to add a "_bcqt" postfix and save in the same directory.
-h, --help|Show this help message and exit.


## Process randomer for Picard on marking PCR duplicates
### 1. Extract randomer from fastq and align to reference
Command example:

`python yl_fastq_tools.py -f input.fastq -s 31:50 -r 1:9 -o $DIR/input_randomer.fastq`

Randomer sequence and corresponding sequencing quality score will be extracted and stored in the sequence identifier (the first line of a read). For example, the following read:

    @M02357:256:000000000-ARWLK:1:1101:23366:8048 1:N:0:1
    CGATGTGCTTGTGGAAAGGACGAAACACCGCTGATGGAATAGGAAGCCGTGTTTAAGAGCTATGCTGGAAACAGCA
    \+
    A1AAA1F3DFF1AAFGEAGFFEEEFAHGEGE0A/BA10FBAA1FFCAGAGGFFFFEAGBGEHFFDGDCGF@FF>GH

will be transformed to:

    @CGATGTGCT:A1AAA1F3D:M02357:256:000000000-ARWLK:1:1101:23366:8048 1:N:0:1
    CTGATGGAATAGGAAGCCGT
    \+
    E0A/BA10FBAA1FFCAGAG

The output fastq can now be aligned using your aligner of choice. For example,

`bowtie2 -L 20 -N 0 --no-1mm-upfront -x reference -U $DIR/input_randomer.fastq -S $DIR/input.sam`

### 2. Attach randomer to the BC&QT tag in the sam file
### 3. Picard MarkDuplicates with BARCODE_TAG option
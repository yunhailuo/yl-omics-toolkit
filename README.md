# yl-omics-toolkit
A collection of bioinformatic tools customized to my needs.

**Table of Contents**
- [Usage](#user-content-usage)
  - [yl_fastq_tools](#user-content-yl_fastq_tools)
  - [yl_sam_tools](#user-content-yl_sam_tools)
- [Process randomer for Picard on marking PCR duplicates](#user-content-process-randomer-for-picard-on-marking-pcr-duplicates)
  - [1. Extract randomer from fastq](#user-content-1-extract-randomer-from-fastq)
  - [2. Attach randomer to the BC&QT tag in the sam file](#user-content-2-attach-randomer-to-the-bcqt-tag-in-the-sam-file)
  - [3. Picard MarkDuplicates with BARCODE_TAG option](#user-content-3-picard-markduplicates-with-barcode_tag-option)

## Usage
### yl_fastq_tools
`yl_fastq_tools.py [-h] -i INPUT [INPUT ...] (-x INDEX | -r RANDOMER) [-o OUTPATH]`
### Arguments
&nbsp;&nbsp;&nbsp;Argument&nbsp;&nbsp;&nbsp;|Help message
---|---
-i|Illumina fastq sequencing file(s). Support multi-file input. A separated result folder and report will be generated for each fastq input.
-x, --index|Basepair position for library indices. Fastq file will be demultiplex accordingly
-r, --randomer|Basepair position for randomer, which indicates unique DNA/RNA molecules before PCR. Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b". For example: "1:9,31:50" stands for a sequence combining nucleotides from position 1 to position 9 with nucleotides from position 31 to position 50.
-o, --outpath|Path for outputs. Default is the same directory as input fastq(s).
-h, --help|Show this help message and exit.
### yl_sam_tools
`yl_sam_tools.py [-h] -s SAM [-o OUTSAM]`
### Arguments
&nbsp;&nbsp;&nbsp;Argument&nbsp;&nbsp;&nbsp;|Help message
---|---
-s --sam|One aligned sam file which needs to attach barcodes at BC&QT tags. The barcode and its sequencing quality should have been stored in the read's name, separated by ":". Check "yl_fastq_tools" module for detailed formats.
-o, --outsam|Output unsorted sam file. Default is to add a "_bcs" postfix and save in the same directory.
-h, --help|Show this help message and exit.


## Process randomer for Picard on marking PCR duplicates
### 1. Extract randomer from fastq
### 2. Attach randomer to the BC&QT tag in the sam file
### 3. Picard MarkDuplicates with BARCODE_TAG option
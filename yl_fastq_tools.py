#!/usr/bin/env python
#
# Copyright (c) 2016 Yunhai Luo
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies ofthe Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""The module contains a few tools useful for handling NGS data.
    The initial aim of this module is to extract molecular barcodes (randomers)
    from fastq sequences and use these info after alignment to remove PCR
    duplicates through BT&QT tags in the sam file and picard's MarkDuplicates
    functionality combinatorially.
"""
import argparse
import itertools


def fastq_reads(fastq):
    """Get reads from fastq file.
    A generator function which iterates over and returns fastq reads by reading
    in 4 lines a time and return the read as a list of 4. Fastq file won't be
    validated. Thus errors won't be detected or handled properly.
    Arg:
        fastq: Path, including name, of one fastq file.
    Yield:
        A list with 4 elements having the expected order of seqence ID,
        sequence, Phred Quality Scores ID, and Phred Quality Scores. For
        example:
        ['@M02357:256:000000000-ARWLK:1:1101:17612:1950 1:N:0:1',
         'TTAGTTTCGTGTGGAAAGGACGAAACACCGCCTCGAGACTCGTTTACTAGGTTTAAGAGCTATGCT',
         '+',
         '>>11>DDD1C1>11111A11B0000A0FE0AEEA////A11/BEF/1A210DDGB2111BF1F1FG']
    Notes:
        Fastq file won't be validated. Therefore, blank or non-fastq line(s) in
        the file may knock off (shift) the expect frame of fastq reads. Also, if
        there are more than one blank lines at the end of file, all but the last
        blank line will be returned as ''.
    """

    with open(fastq) as fastq_file:
        while True:
            read = [line.rstrip('\n') for line in itertools.islice(fastq_file,
                                                                   4)]
            if not read:
                break
            else:
                yield read


def multi_idcs(idcs_str):
    """Convert a string input to indices used by the multi_slices function.
    The main purpose for this function is to generate a ready to use indices
    for the multi_slices function. It saves a little bit time but may save
    quite a bit if the multi_slices function will be used multiple times in a
    loop. It might be better if these two functions are packed in a new class,
    such as "MultiSlicer".
    Arg:
        idcs_str: A string represents the position of substring(s) of interest.
            The format is as follow:
            - One single number stands for a single character of interest;
            - "i:j" stands for a range of interst, staring from position i and
              end at postion j. "i:", ":j", and ":" are all allowed. Check
              "Python - Sequence Types" document for details.
            - Comma is used to separate multiple discrete postions of interest.
              For example: "1:9,31:50".
        This argument will not be validated. Characters other than number, ":",
        and "," will be converted to integer and may lead to unexpected results.
    Return:
        A list of tuples where each tuple represents for a standard python
        sequence slicing operations. For example:
        [(0,9), (30,50)]
        It stands for a string from index 0 to index 8 together with a string
        from index 30 to index 49 (all zero based index). Check the multi_slices
        function for other details.
    """

    idcs = []
    for idx_str in idcs_str.split(','):
        if ':' in idx_str:
            i, j = idx_str.split(':')
            if i == '':
                i = 1
            if j == '':
                idcs.append((int(i)-1, None))
            else:
                idcs.append((int(i)-1, int(j)))
        else:
            idcs.append((int(idx_str)-1, int(idx_str)))
    return idcs


def multi_slices(seq, idcs):
    """Slice out multiple substrings and concatenate them as a new single string
    Args:
        seq: The parent string to be sliced.
        idcs: A list of tuples where each tuple is a pair of indices denoting
            limits of one substring. This can be generated from the multi_idcs
            function. Check it for details.
    Return:
        A string which concatenates all substrings marked by the "idcs" in the
        same order as that in the "idcs".
    """

    sublist = [seq[i:j] for (i, j) in idcs]
    return ''.join(sublist)


def _seqmatch(query, ref, query_allow_n=False):
    """Check if DNA sequences in "query" and "ref" match
    This function uses a slow algorithm and should only be used for a sequencing
    index matching during demultiplexing. "N" in a "ref" sequence will be
    considered as matching to any bases in a "query" sequence, while "N" in a
    "query" sequence will by default be considered as NOT matching to any bases
    in a "ref" sequence.
    Args:
        query: A string of the index sequencing result.
        ref: A string of expected index sequence.
        query_allow_n: By default, a "N" base in a "query" sequence is NOT
            allowed to match any bases in a "ref" sequence. This behavior can be
            changed by switching "query_allow_n" argument to "True".
    Return:
        A Boolean value standing for match or not match.
    """

    # "==" or "!=" should be faster than the following code for exact match.
    # So do it first!
    if query != ref:
        if len(query) != len(ref):
            return False
        for q_char, r_char in zip(query, ref):
            if q_char == r_char:
                continue
            elif r_char == 'N':
                continue
            elif (q_char == 'N') and query_allow_n:
                continue
            else:
                return False
    return True


def _demux_dict(output):
    """Convert a string to a index to file dictionary used by the "fastq_write"
        function during demultiplexing.
    Arg:
        output: A string describing the index to file dictionary. Each
        index-file set will be separated by ";". Within a index-file set, index
        and file name are separated by ":". For index, more than one indices can
        be assigned to the same file. In that case, different indices are
        separated by ",". For example:
        "ATCACG,CGATGT,TTAGGC:sample1.fastq;TGACCA:sample2.fastq;
        ACAGTG,GCCAAT:sample3.fastq"
        For this "output" argument:
            - Reads with ATCACG or CGATGT or TTAGGC index sequence will be
              written to "sample1.fastq";
            - Reads with TGACCA index sequence will be writen to
              "sample2.fastq";
            - Reads with either ACAGTG or GCCAAT will be written to
              "sample3.fastq".
            - Reads having none of listed index sequences will be written to a
              "nonindex_file". Check out the "_fastq_write" function for how to
              specify a "nonindex_file".
    Return:
        A dictionary whose key is an expected index sequence and corresponding
        value is a (fastq) file object to which reads with expected index should
        be written to.
    """

    idx_file_dict = {}
    for file_idcs in output.split(';'):
        idcs, filename = file_idcs.split(':')
        file = open(filename, 'w')
        for idx in idcs.split(','):
            idx_file_dict[idx] = file
    return idx_file_dict


def _fastq_write(read, nonindex_file, read_index=None, index_file_dict=None):
    """Write reads into target fastq file(s)
    The main purpose for this function is to demultiplex fastq file according to
    sequencingn index. Though it's not recommended, it can be
    configured/simplified to write a read into one fastq file without any index
    references (give only "read" and "nonindex_file" arguments; use
    "nonindex_file" object as the output fastq).
    Args:
        read: A list of 4 strings in the order of seqence ID, sequence, Phred
            Quality Scores ID, and Phred Quality Scores. This read will be write
            into a fastq file as it is; in other words, this read will not be
            proceesed.
        read_index: A string of the sequencing index for the input "read". It
            will be used to identify the corresponding output fastq file
            (demultiplex).
        index_file_dict: A dict generated by the "_demux_dict" function. Reads
            have no expected index will be written to the fastq output specified
            by "nonindex_file". Check out the "_demux_dict" function for
            details.
        nonindex_file: A file object to which any reads whose index fail to
            match the "index_file_dict" will be written. It can also be used to
            perform a simple fastq file writing without demultiplexing. In fact,
            if a "read_index" is not provided, "index_file_dict" will be ignored
            and read will only be written to this "nonindex_file".
    """

    # Default output is the nonindex_file. Demultiplex when both read_index and
    # index_file_dict are provided; otherwiseIf write read to the nonindex_file.
    output_fastq = nonindex_file
    if read_index and index_file_dict:
        for expect_index in index_file_dict:
            # The following sequence match function allow "N"
            if _seqmatch(read_index, expect_index):
                output_fastq = index_file_dict[expect_index]
                break
    read.append('') # Add the last line break
    output_fastq.write('\n'.join(read))


def process_fastq(fastq, bc_bp=None, slice_bp=None, index_bp=None, output=None):
    """Process fastq file
    For each read, this function can clip sequence, extracting molecular
    barcodes (randomer) and prepend them to the sequence id, and split it to
    individual fastq output (demultiplex).
    Args:
        fastq: Input with raw data from Illuminar next-gen sequencing.
        bc_bp: Base pair position for molecular barcodes (randomer). It is
            1-based and allows discrete regions. Check the "multi_idcs" function
            for details. Example: "1:9,31:50".
        slice_bp: Base pair position for sequence of interest to keep in output.
            It is 1-based and allows discrete regions. Check the "multi_idcs"
            function for details. Example: "1:9,31:50".
        index_bp: Base pair position for sequencing index, if a custom
            demultiplexing is needed.
        output: If no demultiplexing is need, this argument will accept a path
        for one single fastq output. If not provided, the default output file
        will be saved in the same directory as the fastq input with a  "_noidx"
        postfix. If a custom demultiplexing is needed (with index_bp argument
        provided), this argument will accept a string describing the index to
        file dictionary. Each index-file set will be separated by ";". Within a
        index-file set, index and file name are separated by ":". For index,
        more than one indices can be assigned to the same file. In that case,
        different indices are separated by ",". For example:
        "ATCACG,CGATGT,TTAGGC:sample1.fastq;TGACCA:sample2.fastq;
        ACAGTG,GCCAAT:sample3.fastq"
        For this "output" argument:
            - Reads with ATCACG or CGATGT or TTAGGC index sequence will be
              written to "sample1.fastq";
            - Reads with TGACCA index sequence will be writen to
              "sample2.fastq";
            - Reads with either ACAGTG or GCCAAT will be written to
              "sample3.fastq".
            - Reads having none of listed index sequences will be written to the
              a default output file which will be saved in the same directory as
              the fastq input with a  "_noidx" postfix.
    """

    # Prepare indices and output file(s)
    if bc_bp:
        bc_idcs = multi_idcs(bc_bp)
    if slice_bp:
        slice_idcs = multi_idcs(slice_bp)
    final_out = fastq.replace('.fastq', '_noidx.fastq')
    out_index = None
    output_dict = None
    if index_bp:
        index_idcs = multi_idcs(index_bp)
        output_dict = _demux_dict(output)
    elif output:
        final_out = output
    nonindex_file = open(final_out, 'w')

    # Process reads
    for read in fastq_reads(fastq):
        seq_id, seq, _, pqs = read
        if bc_bp:
            read[0] = '@{}'.format(':'.join([multi_slices(seq, bc_idcs),
                                             multi_slices(pqs, bc_idcs),
                                             seq_id[1:]]))
        if slice_bp:
            read[1] = multi_slices(seq, slice_idcs)
            read[3] = multi_slices(pqs, slice_idcs)
        if index_bp:
            out_index = multi_slices(seq, index_idcs)
        _fastq_write(read, nonindex_file, out_index, output_dict)

    # Close all files
    nonindex_file.close()
    if index_bp:
        for file in set(output_dict.values()):
            file.close()


def main():
    """Use the main() function to test the major function of this module:
        - Process fastq file
    """

    parser = argparse.ArgumentParser(description='Pre-procees one fastq file '
                                                 'for downstream analysis.')
    parser.add_argument('-f', '--fastq',
                        help='One fastq file to be processed',
                        required=True)
    parser.add_argument('-s', '--sequence',
                        help='Basepair position for target sequences which '
                             'should be kept for downstream analysis. Support '
                             'both continuous and gapped regions. Use syntax '
                             '\'start1:stop1,start2:stop2,...\'. Both the start'
                             ' and the end are included. If a stop is omitted, '
                             'start-stop pairs after it will be ignored and '
                             'sequencing bases will be kept up to the end. If '
                             'this option is omitted, the whole sequence will '
                             'be kept. However, at least one of the \'-s\', '
                             '\'-i\'+\'-b\' and \'-r\' options must be '
                             'provided.',
                        default=None)
    parser.add_argument('-i', '--barcode-index',
                        help='Basepair position for library barcodes which help'
                             ' demultiplex reads. Barcode sequences and '
                             'corresponding output file names must be provided'
                             ' in the \'-b\' option (check it for details). '
                             'Support both continuous and gapped regions. Use '
                             'syntax \'start1:stop1,start2:stop2,...\'. The '
                             'syntax specification is the same as \'-s\' option'
                             ' above (check it for details). At least one of '
                             'the \'-s\', \'-i\'+\'-b\' and \'-r\' options must'
                             ' be provided.',
                        default=None)
    parser.add_argument('-b', '--barcode-file',
                        help='Barcode sequences and corresponding output file '
                             'names for demultiplexing. Use syntax '
                             '\'barcode1-1,barcode1-2,barcode1-3,...:filename1;'
                             'barcode2-1,barcode2-2,barcode2-3,...:filename2;'
                             '...\'. No space is allowed. "N" is supported in '
                             'barcodes but not in read sequences. At least one '
                             'of the \'-s\', \'-i\'+\'-b\' and \'-r\' options '
                             'must be provided.',
                        default=None)
    parser.add_argument('-r', '--randomer',
                        help='Basepair position for randomer '
                             '(molecular barcode) which help remove PCR '
                             'duplicates generated during library preparation. '
                             'It will be extracted and saved in corresponding '
                             'reads\' name temporarily. Later in the pipeline '
                             '(in the sam file), it can be moved to the BC&QT '
                             'tag which can be utilized by the Picard tool for '
                             'marking/removing duplicates. Support both '
                             'continuous and gapped regions. Use syntax '
                             '\'start1:stop1,start2:stop2,...\'. The syntax '
                             'specifications are the same as \'-s\' option '
                             'above (check it for details). At least one of the'
                             ' \'-s\', \'-i\'+\'-b\' and \'-r\' options must be'
                             ' provided.',
                        default=None)

    args = parser.parse_args()
    if args.subcom == 'process_fastq':
        process_fastq(args.fastq, bc_bp=args.barcode, slice_bp=args.slice,
					  index_bp=args.index, output=args.outfastq)
    else:
        parser.print_help()
        exit(1)


if __name__ == '__main__':
    main()


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
    The initial aim of this module is to attach molecular barcodes (randomers)
    extracted during fastq processing (as part of the sequence name) to the
    BT&QT tag field in the sam file so that the picard tool can MarkDuplicates.
"""
import argparse
import itertools
import pysam


def attach_barcode(sam, output):
    """Attach molecular barcodes (randomers) to the BC&QT tag in the sam file.
    Args:
        sam: A sam file which needs to attach barcodes at BC&QT tags. The
            barcode and its sequencing quality should be stored in the read's
            name, separated by ":". Check the "process_fastq" function for
            details.
        output: Path for a unsorted sam file output. Default path is the current
            working directory (os.getcwd()) and the default name is to add a
            "_bcqt" postfix to the name of the input sam file.
    """
    if output is None:
        output = sam.replace('.sam', '_bcqt.sam')
    infile = pysam.AlignmentFile(sam, "r")
    outfile = pysam.AlignmentFile(output, "wh", template=infile)
    for read in infile.fetch():
        id_sam = read.query_name
        sep_si = id_sam.index(':')
        bc_seq = id_sam[0:sep_si]
        sep_qi = sep_si + 1 + len(bc_seq)
        bc_pqs = id_sam[sep_si + 1: sep_qi]
        read.set_tag('BC', bc_seq)
        read.set_tag('QT', bc_pqs)
        read.query_name = id_sam[sep_qi+1:]
        outfile.write(read)
    outfile.close()
    infile.close()


def main():
    """Use the main() function to test major functions of this toolkit:
        - Process fastq file
        - Attach molecular barcodes (randomers) to the BC&QT tag of a sam file.
    """

    parser = argparse.ArgumentParser(description='Yunhai Luo\'s toolkit')
    subparsers = parser.add_subparsers(dest='subcom')

    # Sub-command for attaching randomer sequences to the BC&QT tag, which will
    # be used by picard for "MarkDuplicates".
    parser_attach = subparsers.add_parser('attach_sam',
                                          help='Attach barcodes to the BC&QT '
                                               'tag in one aligned sam file.')
    parser_attach.add_argument('-s', '--sam',
                               help='One aligned sam file which needs to attach'
                                    ' barcodes at BC&QT tags. The barcode and '
                                    'its sequencing quality should be stored in'
                                    'the read\'s name, separated by ":". Check '
                                    '"process_fastq" sub-command for details.',
                               required=True)
    parser_attach.add_argument('-o', '--outsam',
                               help='Output unsorted sam file. Default is to '
                                    'add a "_bcs" postfix and save in the same'
                                    ' directory.',
                               default=None)

    args = parser.parse_args()
    if args.subcom == 'attach_sam':
        attach_barcode(args.sam, args.outsam)
    else:
        parser.print_help()
        exit(1)


if __name__ == '__main__':
    main()


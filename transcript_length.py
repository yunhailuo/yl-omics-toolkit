#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Yunhai Luo
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to 
# deal in the Software without restriction, including without limitation the 
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
# sell copies ofthe Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
# IN THE SOFTWARE.

"""
This module provides a function for calculating non-overlapping total exon 
length of gene from a GTF/GFF2 or GFF3 annotation file. The generated gene 
length table can be used for calculating RPKM/FPKM/TPM.

This module can also be used as a script with one optional option (``-o`` for 
output table) and a required argument (for input GFF file). Please run 
``python transcript_length.py -h`` for details.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import print_function

import argparse
import os
import re

import pandas as pd

def gff_transcript_length(path, output):
    """Calculate non-overlapping total exon length for each gene from a 
    GTF/GFF2 or GFF3 annotation file.
    
    Args:
        path (str, pathlib.Path, py._path.local.LocalPath or any object with a 
        read() method): One GTF/GFF2 or GFF3 file. This input file will be 
            read by ``pandas.read_table``. Therefore, it can any type accepted 
            by ``pandas.read_table``.
        output (str or filehandle): File path, including both directory and 
        filename, for output transcript length table. It will be passed to 
        ``pandas.DataFrame.to_csv``.
    """
    
    _, ext = os.path.splitext(path)
    print('Loading {} file ...'.format(ext[1:].upper()), end='')
    df = pd.read_table(path, names=['chrom', 'feature', 'chromStart', 
                                    'chromEnd', 'attr'], 
                       usecols=[0, 2, 3, 4, 8], comment='#')
    print('\rFiltering exons ...', end='')
    df = df[df['feature']=='exon']
    print('\rExtracting Ensembl IDs ...', end='')
    if ext.lower() == '.gtf' or ext.lower() == '.gff2':
        gene_id_parser = lambda x: re.search(
                r'gene_id "([^"]+)"(;|$)', x
            ).group(1)
    elif ext.lower() =='.gff3':
        gene_id_parser = lambda x: re.search(
                r'gene_id=([^;\n]+)(;|$)', x
            ).group(1)
    else:
        raise TypeError('Unknown file extension: {}'.format(ext))
    df['gene_id'] = df['attr'].apply(gene_id_parser)
    print('\rSorting the table for transcript length calculation ...', end='')
    df = df.drop(['feature', 'attr'], axis=1).sort_values(
            ['gene_id', 'chrom', 'chromStart', 'chromEnd']
        )
    print('\rCalculating non-overlapping exon length for each gene ...', 
          end='')
    (df.groupby(['gene_id', 'chrom'])
       .apply(nonoverlap_len)
       .reset_index()
       .drop('chrom', axis=1)
       .to_csv(output, sep='\t', header=False, index=False))
    print('\rTranscript length table ready.                           ')


def nonoverlap_len(bed):
    """Calculate total non-overlap length of a set of regions.
    
    Args:
        bed (pandas.DataFrame): A sorted dataframe in BED format. It should 
            have at least four columns for "chrom", "chromStart", "chromEnd" 
            and "geneName/ID". This sorted dataframe must have a single 
            "chrom" value for all regions and is expected to be sorted 
            first by "chromStart" and then by "chromEnd" Columns should be in 
            the expect order (i.e. values will be selected by position, 
            ``iloc``), though the name is not important.
    
    Return:
        length (int): the total non-overlap length of all regions in the input 
            BED.
    """
    
    assert len(bed.iloc[:, 0].unique()) == 1
    length = 0
    start = 0
    end = -1
    for _, row in bed.iterrows():
        if row.iloc[1] > end + 1:
            length += end - start + 1
            start = row.iloc[1]
            end = row.iloc[2]
        elif row.iloc[2] > end:
            end = row.iloc[2]
    length += end - start + 1
    return length


def main():
    parser = argparse.ArgumentParser(
            description='Get a table of gene length for all transcripts from '
                        'a GTF/GFF2/GFF3 file, which can be potentially used '
                        'for calculating RPKM/FPKM.'
        )
    parser.add_argument('gff', metavar='GFF', type=str,
                        help='One input GTF/GFF2/GFF3 file from which the '
                             'length of transcript for each gene will be '
                             'calculated.')
    parser.add_argument('-o', '--output', help='Optional. Specify path and '
                        'filename for output gene length table. Default '
                        'saving directory is that of the input file and '
                        'default filename is "transcript_length.tsv".',
                        default='transcript_length.tsv')
    args = parser.parse_args()
    os.chdir(os.path.dirname(os.path.abspath(args.gff)))
    gff_transcript_length(args.gff, args.output)


if __name__ == '__main__':
    main()

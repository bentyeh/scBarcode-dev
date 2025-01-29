import argparse
import collections
import os
import re
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import positive_int, grouper

import pandas as pd
import pysam


def main():
    args = parse_arguments()
    dedup_fun = dedup_paired_end if args.paired else dedup_single_end
    dedup_fun(
        args.input,
        path_out_bam=args.output,
        path_out_bed=args.counts,
        barcode_rgx=args.barcode_rgx,
        threads=args.threads
    )


def dedup_single_end(
    path_in_bam: str,
    path_out_bam: str | None = None,
    path_out_bed: str | None = None,
    barcode_rgx: str | None = None,
    threads: int = 1
) -> pd.DataFrame:
    '''
    Deduplicate aligned, single-end reads by genomic coordinates and optional barcode.
    - Unlike samtools markdup, this function uses coordinates of both the leftmost and
      rightmost mapped bases.

    Args
    - path_in_bam: path to name-sorted BAM file
    - path_out_bam: path to output BAM file of deduplicated reads. If None, write to standard out.
    - path_out_bed: path to output BED file of read counts.
    - barcode_rgx: Regular expression for barcode in the read name
        Currently only supports 1 capture group for an integer.
    - threads: Number of threads to use for reading and writing BAM files

    Returns: Pandas DataFrame of read counts
        Columns = chr, start, end, barcode, count
        Coordinates are 0-based (BED format).
    '''
    if barcode_rgx:
        regex_barcode = re.compile(barcode_rgx)
    path_out_bam = path_out_bam if path_out_bam is not None else sys.stdout.buffer
    path_in_bam = path_in_bam if path_in_bam != '-' else sys.stdin.buffer

    entries = collections.Counter()
    with pysam.AlignmentFile(path_in_bam, 'rb', threads=threads) as file_in:
        header = file_in.header.to_dict()
        with pysam.AlignmentFile(path_out_bam, 'wb', threads=threads, header=header) as file_out:
            for read in file_in.fetch(until_eof=True):
                barcode = int(regex_barcode.search(read.qname).groups()[0]) if barcode_rgx else '-'
                entry = (read.reference_id, read.reference_start, read.reference_end, barcode)
                if entry not in entries:
                    file_out.write(read)
                entries[entry] += 1
    s = pd.Series(entries, name='count')
    s.index.set_names(['chr', 'start', 'end', 'barcode'], inplace=True)
    df = s.reset_index()
    REFID_TO_CHR = {i: SQ['SN'] for i, SQ in enumerate(header['SQ'])}
    DTYPE_REFID = pd.CategoricalDtype(categories=list(range(len(header['SQ']))), ordered=True)
    df['chr'] = df['chr'].astype(DTYPE_REFID).cat.rename_categories(REFID_TO_CHR)
    df.sort_values(['chr', 'start', 'end', 'barcode', 'count'], inplace=True)
    if path_out_bed:
        df.to_csv(path_out_bed, sep='\t', index=False, header=False)
    return df


def dedup_paired_end(
    path_in_bam: str,
    path_out_bam: str | None,
    path_out_bed: str | None,
    barcode_rgx: str | None = None,
    threads: int = 1
):
    '''
    Deduplicate aligned, paired-end reads by genomic coordinates and optional barcode.
    - Unlike samtools markdup, read orientation is not considered.

    Args
    - path_in_bam: path to name-sorted BAM file
    - path_out_bam: path to output BAM file of deduplicated reads. If None, write to standard out.
    - path_out_bed: path to output sorted BED file of read counts.
    - barcode_rgx: Regular expression for barcode in the read name
        Currently only supports 1 capture group for an integer.
    - threads: Number of threads to use for reading and writing BAM files

    Returns: Pandas DataFrame of read counts
        Columns = chr, start, end, barcode, count
        Coordinates are 0-based (BED format).
    '''
    if barcode_rgx:
        regex_barcode = re.compile(barcode_rgx)
    path_out_bam = path_out_bam if path_out_bam is not None else sys.stdout.buffer
    path_in_bam = path_in_bam if path_in_bam != '-' else sys.stdin.buffer

    entries = collections.Counter()
    with pysam.AlignmentFile(path_in_bam, 'rb', threads=threads) as file_in:
        header = file_in.header.to_dict()
        with pysam.AlignmentFile(path_out_bam, 'wb', threads=threads, header=header) as file_out:
            for read1, read2 in grouper(file_in.fetch(until_eof=True), 2, incomplete='strict'):
                assert read1.qname == read2.qname
                assert read1.reference_name == read2.reference_name
                assert read1.reference_end >= read1.reference_start
                assert read2.reference_end >= read2.reference_start
                assert read1.template_length == -read2.template_length
    
                barcode = int(regex_barcode.search(read1.qname).groups()[0]) if barcode_rgx else '-'
                
                if read1.is_reverse:
                    assert read2.is_forward
                    entry = (read1.reference_id, read2.reference_start, read1.reference_end, barcode)
                else:
                    assert read1.is_forward and read2.is_reverse
                    entry = (read1.reference_id, read1.reference_start, read2.reference_end, barcode)
                assert entry[2] >= entry[1]
                assert entry[2] - entry[1] == abs(read1.template_length)
                if entry not in entries:
                    file_out.write(read1)
                    file_out.write(read2)
                entries[entry] += 1
    s = pd.Series(entries, name='count')
    s.index.set_names(['chr', 'start', 'end', 'barcode'], inplace=True)
    df = s.reset_index()
    REFID_TO_CHR = {i: SQ['SN'] for i, SQ in enumerate(header['SQ'])}
    DTYPE_REFID = pd.CategoricalDtype(categories=list(range(len(header['SQ']))), ordered=True)
    df['chr'] = df['chr'].astype(DTYPE_REFID).cat.rename_categories(REFID_TO_CHR)
    df.sort_values(['chr', 'start', 'end', 'barcode', 'count'], inplace=True)
    if path_out_bed:
        df.to_csv(path_out_bed, sep='\t', index=False, header=False)
    return df


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Remove duplicate reads based on identical genomic alignment coordinates."
    )
    parser.add_argument(
        "input",
        metavar="in.bam|-",
        help="Input BAM file. Use '-' for standard in."
    )
    parser.add_argument(
        "-o", "--output",
        metavar="out.bam",
        help="Output BAM file. If not provided, write to standard out."
    )
    parser.add_argument(
        "-c", "--counts",
        metavar="counts.bed(.gz)",
        help="Output counts BED file. Columns = chr, start, end, barcode, count."
    )
    parser.add_argument(
        "-p", "--paired",
        action="store_true",
        help="Input is a name-sorted paired-end alignment file."
    )
    parser.add_argument(
        "-t", "--threads",
        type=positive_int,
        default=1,
        metavar="#",
        help="Number of threads to use for compressing/decompressing BAM files",
    )
    parser.add_argument(
        "--barcode-rgx",
        metavar="REGEX",
        help=("Regular expression for barcode in the read name. Identify duplicates by "
              "alignment coordinates and barcode.")
    )
    return parser.parse_args()


if __name__ == '__main__':
    main()
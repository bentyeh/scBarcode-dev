import argparse
import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from helpers import positive_int

import pysam


def main():
    args = parse_arguments()
    remove_unpaired(
        args.input,
        path_out_bam=args.output,
        threads=args.threads
    )


def remove_unpaired(
    path_in_bam: str,
    path_out_bam: str | None = None,
    threads: int = 1
) -> None:
    '''
    Remove unpaired reads by read name.

    Args
    - path_in_bam: path to name-collated BAM file
    - path_out_bam: path to output BAM file of deduplicated reads. If None, write to standard out.
    - threads: Number of threads to use for reading and writing BAM files

    Returns: None
    '''
    path_out_bam = path_out_bam if path_out_bam is not None else sys.stdout.buffer
    path_in_bam = path_in_bam if path_in_bam != '-' else sys.stdin.buffer

    in_pair = False
    paired_read = None
    with pysam.AlignmentFile(path_in_bam, 'rb', threads=threads) as file_in:
        header = file_in.header.to_dict()
        with pysam.AlignmentFile(path_out_bam, 'wb', threads=threads, header=header) as file_out:
            for read in file_in.fetch(until_eof=True):
                if in_pair and read.qname == paired_read.qname:
                    file_out.write(paired_read)
                    file_out.write(read)
                    in_pair = False
                else:
                    paired_read = read
                    in_pair = True


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Remove unpaired reads based on identical genomic alignment coordinates."
    )
    parser.add_argument(
        "input",
        metavar="in.bam|-",
        help=("Input BAM file, with reads grouped by name (i.e., after running "
              "samtools collate) such that reads from a read pair are adjacent. "
              "Use '-' for standard in.")
    )
    parser.add_argument(
        "-o", "--output",
        metavar="out.bam",
        help="Output BAM file. If not provided, write to standard out."
    )
    parser.add_argument(
        "-t", "--threads",
        type=positive_int,
        default=1,
        metavar="#",
        help="Number of threads to use for compressing/decompressing BAM files",
    )
    return parser.parse_args()


if __name__ == '__main__':
    main()
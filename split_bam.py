#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input bam file")
parser.add_argument("--prefix", default="chunk")
parser.add_argument("--chunk_size", type=int, default=50000)
parser.add_argument("-d", action="store_true", default=False)
args = parser.parse_args()

import pysam

infile = pysam.AlignmentFile(args.infile)

chunk_size = args.chunk_size
outfile_pattern = args.prefix + ".{}.bam"

chunk = 0
reads_in_this_chunk = 0
outfile = pysam.AlignmentFile(outfile_pattern.format(chunk), "wb", template=infile)

for read in infile.fetch(until_eof=True):

    # write to output
    outfile.write(read)
    reads_in_this_chunk += 1

    # print progress
    if reads_in_this_chunk % 100 == 0:
        sys.stderr.write(
            "\r{} written to {}".format(
                reads_in_this_chunk, outfile_pattern.format(chunk)
            )
        )

    # switch output file
    if reads_in_this_chunk >= chunk_size:
        outfile.close()
        reads_in_this_chunk = 0
        chunk += 1
        outfile = pysam.AlignmentFile(
            outfile_pattern.format(chunk), "wb", template=infile
        )
        sys.stderr.write("\n")


outfile.close()

#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument(
    "infile",
    nargs="?",
    help="input bam file",
    type=argparse.FileType("r"),
    default=sys.stdin,
)
parser.add_argument(
    "outfile",
    nargs="?",
    help="output bam file",
    type=argparse.FileType("w"),
    default=sys.stdout,
)
parser.add_argument("--minlength", "-l", type=int, default=6000)
args = parser.parse_args()

import pysam


sam = pysam.AlignmentFile(args.infile)
out = pysam.AlignmentFile(args.outfile, "wb", template=sam)

for read in sam.fetch(until_eof=True):
    if read.reference_length > args.minlength:
        out.write(read)

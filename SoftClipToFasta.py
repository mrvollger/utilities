#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", help="positional input")
parser.add_argument(
    "-d", help="store args.d as true if -d", action="store_true", default=False
)
args = parser.parse_args()

import pysam


# soft clip
S = 4

bam = pysam.AlignmentFile(args.infile)

for rec in bam.fetch(until_eof=True):
    cigar_ends = (rec.cigartuples[0], rec.cigartuples[-1])
    for idx, (opt, length) in enumerate(cigar_ends):
        if opt == S:
            seq = ""
            if idx == 0:
                seq = rec.seq[0:length]
            elif idx == 1:
                seq = rec.seq[-length:]

            print(
                ">{}.{} soft clipped length {}\n{}".format(
                    rec.query_name, idx + 1, length, seq
                )
            )

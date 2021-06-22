#!/usr/bin/env python
import argparse
import os
import sys
import pysam

# global var for inputs
args = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="positional input")
    parser.add_argument(
        "-n", help="length of each window", type=int, default=5 * 10 ** 6
    )
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    args = parser.parse_args()

    fasta = pysam.FastaFile(args.infile)

    for name, length in zip(fasta.references, fasta.lengths):
        seq = fasta.fetch(name)

        for start in range(0, length, args.n):
            end = min(length, start + args.n)
            sys.stdout.write(
                ">{}:{}-{}\n{}\n".format(name, start + 1, end, seq[start:end])
            )

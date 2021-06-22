#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import os
import sys
import argparse
import pysam


def find_gc(name, seq):
    length = len(seq)
    i = 1
    while i < len(seq):
        pre = seq[i - 1]
        cur = seq[i]
        if pre == "C" and cur == "G":
            print(name + "\t" + str(i - 1) + "\t" + str(i + 1) + "\t+")
        if pre == "G" and cur == "C":
            print(name + "\t" + str(i - 1) + "\t" + str(i + 1) + "\t-")
        i += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="positional input")
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    args = parser.parse_args()
    fasta = pysam.FastaFile(args.infile)
    for name in fasta.references:
        seq = fasta.fetch(name)
        find_gc(name, seq)

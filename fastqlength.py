#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("fastq", nargs="+", help="input fastq file")
args = parser.parse_args()

import gzip
from Bio import SeqIO
import sys

lengths = {}


def getlens(handle):
    # print(handle, file=sys.stderr)
    for record in SeqIO.parse(handle, "fastq"):
        tmp = len(record.seq)
        if tmp not in lengths:
            lengths[tmp] = 0
        lengths[tmp] += 1


for fastq in args.fastq:
    if fastq[-3:] == ".gz":
        with gzip.open(fastq, "rt") as handle:
            getlens(handle)
    else:
        with open(fastq, "r") as handle:
            getlens(handle)

for key in sorted(lengths.keys()):
    print("{}\t{}".format(key, lengths[key]))

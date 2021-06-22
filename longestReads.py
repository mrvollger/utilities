#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("bam", help="bam file with the reads")
parser.add_argument(
    "-o",
    "--outbam",
    help="name of outbam, defulats to topLength.bam",
    default="topLength.bam",
)
parser.add_argument(
    "--percent", "-p", help="percent of reads to grab", type=float, default=0.1
)
parser.add_argument("--contig", help="to look at", default=None)
parser.add_argument("-d", action="store_true", default=False)
args = parser.parse_args()
DEBUG = args.d

import glob
import os
import sys
import re
import itertools
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam
import runCmd

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

if not os.path.exists(args.bam + "bai"):
    pysam.index(args.bam, args.bam + ".bai")
bamfile = pysam.AlignmentFile(args.bam, "rb")
topLen = pysam.AlignmentFile(args.outbam, "wb", template=bamfile)


lengths = []
for read in bamfile.fetch():
    lengths.append(read.infer_query_length())

lengths = np.sort(np.array(lengths))
# print a histogram of lengths
length = len(lengths)
for i in range(1, 11):
    posStart = (i - 1) * length / 10
    posEnd = i * length / 10 - 1
    print(
        "{}%-{}%: {}-{}".format(
            (i - 1) * 10, i * 10, lengths[posStart], lengths[posEnd]
        )
    )

minLength = lengths[int((1.0 - args.percent) * length - 1)]
print("Minimum Length = {}".format(minLength))

for read in bamfile.fetch():
    tempLen = read.infer_query_length()
    if tempLen >= minLength:
        topLen.write(read)


topLen.close()
bamfile.close()

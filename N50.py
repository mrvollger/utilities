#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(
    description="This is similar to merge in bedtools, but instead of mergeing overlapping regions, it collapses all regions that are entirly contained within another region"
)
parser.add_argument("infile", nargs="?", type=argparse.FileType("r"), default=sys.stdin)
parser.add_argument(
    "outfile", nargs="?", type=argparse.FileType("w"), default=sys.stdout
)
parser.add_argument("-d", action="store_true", default=False)
parser.add_argument("-n", "--colNum", help="col with read lengths", type=int, default=0)
args = parser.parse_args()

import pandas as pd

table = pd.read_csv(args.infile, header=None, delimiter=r"\s+")
table.sort_values(by=args.colNum, ascending=False, inplace=True)

totalSize = table[args.colNum].sum()


soFar = 0
for readLength in table[args.colNum]:
    soFar += readLength
    if soFar >= totalSize / 2.0:
        print(readLength)
        break

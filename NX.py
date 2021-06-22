#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("infile", nargs="+", help="fai file to get n50 from")
parser.add_argument("--outfile", type=argparse.FileType("w"), default=sys.stdout)
parser.add_argument(
    "-g",
    help="calculate NG50, provide genome size, for hg use 3098794149 ",
    type=int,
    default=None,
)
parser.add_argument("-d", action="store_true", default=False)
parser.add_argument(
    "-c", "--colNum", help="column with readlength in it", type=int, default=2
)
args = parser.parse_args()

import pandas as pd

args.colNum = args.colNum - 1


for infile in args.infile:
    table = pd.read_csv(infile, header=None, delimiter=r"\s+")
    table.sort_values(by=args.colNum, ascending=False, inplace=True)
    totalSize = table[args.colNum].sum()
    if args.g is not None:
        totalSize = args.g

    soFar = 0
    N50 = 0
    for readLength in table[args.colNum]:
        soFar += readLength
        if soFar >= totalSize / 2.0:
            N50 = readLength
            break

    print(infile)
    print(
        "N50_Mbp\tTotal_Gbp\tN50\tTotal_bp\n{:.4f}\t{:.4f}\t{}\t{}".format(
            N50 / 10 ** 6, totalSize / 10 ** 9, N50, totalSize
        )
    )
    if args.d:
        pd.set_option("display.float_format", lambda x: "%.2f" % x)
        print(table[args.colNum].describe())

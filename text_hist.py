#!/usr/bin/env python
import argparse
import os
import sys
import numpy as np

# global var for inputs
args = None


def read_data(f):
    data = []
    # if(f=="-"):
    #    f=sys.stdin
    # else:
    #    f = open(f)
    for idx, line in enumerate(f):
        if idx < args.skip:
            continue
        if line[0] == args.comment:
            continue
        data.append(float(line.strip().split()[args.column - 1]))
    return data


def make_hist(data, args):
    data = np.array(data)

    # get maxes in regular and log space
    MAX = data.max()
    MIN = data.min()

    bins = np.linspace(MIN, MAX, args.bins)
    if args.log:
        LMIN = int(np.floor(np.log10(MIN)))
        LMAX = int(np.ceil(np.log10(MAX)))
        # print(LMIN, LMAX)
        bins = np.logspace(LMIN, LMAX, LMAX - LMIN + 1)
    if args.binwidth:
        bins = np.arange(np.floor(MIN), np.ceil(MAX) + args.binwidth, args.binwidth)
        #print(bins)

    inds = np.digitize(data, bins)
    maxbincount = np.max(np.bincount(inds))

    # print(inds, inds.shape, data.shape, bins.shape, inds.min(), inds.max(), maxbincount)
    fmt = "[{:<15}{:<15}){:>10}\t{}\n"
    sys.stdout.write(
        fmt.format(
            "start", "end", "count", "histogram: *=" + str(maxbincount / args.width)
        )
    )

    for b in range(bins.shape[0]):
        count = (inds == b + 1).sum()
        n = count / maxbincount * args.width
        text = "*" * int(n)
        if int(n) - n != 0:
            text += "."
        #print(b, n, count)
        if b + 1 >= bins.shape[0]:
            sys.stdout.write(
                fmt.format(round(bins[b], 2), round(bins[b], 2), count, text)
            )
        else:
            sys.stdout.write(
                fmt.format(round(bins[b], 2), round(bins[b+1], 2), count, text)
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # parser.add_argument("infiles", nargs="+", help="positional input")
    parser.add_argument("infile", type=argparse.FileType("r"))
    parser.add_argument(
        "-c", "--column", help="column to read from", type=int, default=1
    )
    parser.add_argument("-s", "--skip", help="skip first x lines", type=int, default=0)
    parser.add_argument(
        "--comment", help="skip lines starting with char", type=str, default="#"
    )
    parser.add_argument(
        "-b", "--bins", help="number of bins to make", type=int, default=30
    )
    parser.add_argument(
        "-w", "--binwidth", help="width of bins", type=int, default=None
    )
    parser.add_argument("--width", help="number of bins to make", type=int, default=80)
    parser.add_argument(
        "-l", "--log", help="use log10 base bins", action="store_true", default=False
    )
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    args = parser.parse_args()

    # nargs = len(args.infiles)

    for idx, h in enumerate([args.infile]):
        data = read_data(h)
        make_hist(data, args)

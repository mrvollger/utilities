#!/usr/bin/env python
import argparse
import os
import sys
from numba import njit
from numba import jit
import numpy as np
from scipy.ndimage.measurements import label
import pandas as pd
import networkx as nx


@njit
def intersect(row1, row2, dist=0, limit=2):
    a, a0, a1 = row1
    b, b0, b1 = row2
    if limit:
        dist = min(dist, limit * (a1 - a0 + b1 - b0))
    return (a == b) and ((a1 + dist) >= b0) and ((b1 + dist) >= a0)


@njit
def combos(rows, dist):
    pairs = []
    for i in range(rows.shape[0]):
        row1 = rows[i]
        for j in range(i, rows.shape[0]):
            row2 = rows[j]
            # intersect first range
            if intersect(row1[0:3], row2[0:3], dist=dist):
                # intersect second range
                if intersect(row1[3:6], row2[3:6], dist=dist):
                    pairs.append((i, j))
            else:
                break
    return pairs


@njit
def merge(rows):
    """
    Get the merged row from a set of overlaping bed pairs
    """
    row = np.array(
        (
            rows[0, 0],
            rows[:, 1].min(),
            rows[:, 2].max(),
            rows[0, 3],
            rows[:, 4].min(),
            rows[:, 5].max(),
        )
    )
    return row


def new_rows(rows, pairs):
    g = nx.Graph()
    g.add_edges_from(pairs)
    rtn = []
    for group in nx.connected_components(g):
        group = np.array(list(group))
        to_merge = rows[np.array(group)]
        rtn.append(merge(to_merge))
    return rtn


def write_rows(rows, chrs, out=sys.stdout, symetric=False):
    for r in rows:
        line = (5 * "{}\t" + "{}\n").format(
            chrs[r[0]], r[1], r[2], chrs[r[3]], r[4], r[5]
        )
        out.write(line)
        if symetric:
            line = (5 * "{}\t" + "{}\n").format(
                chrs[r[3]], r[4], r[5], chrs[r[0]], r[1], r[2]
            )
            out.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="positional input")
    parser.add_argument(
        "-d", "--dist", help="distance allowed between intervals", type=int, default=0
    )
    parser.add_argument(
        "-s",
        "--symetric",
        help="report a,b and b,a in the output",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-c",
        "--columns",
        nargs="+",
        help="in order list of columns to use for second set of intervals",
        default=[4, 5, 6],
    )
    # parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
    args = parser.parse_args()
    outfile = sys.stdout
    chr2c, st2c, en2c = [x - 1 for x in args.columns]
    sys.stderr.write("Reading file\n")
    rows = []
    seen = set()
    chr_to_int = {}
    chrs = []
    num = -1
    for line in open(args.infile).readlines():
        if line[0] == "#":
            continue
        t = line.strip().split()
        for c in (t[0], t[3]):
            if c not in seen:
                num += 1
                seen.add(c)
                chrs.append(c)
                chr_to_int[c] = num

        row = (
            chr_to_int[t[0]],
            int(t[1]),
            int(t[2]),
            chr_to_int[t[3]],
            int(t[4]),
            int(t[5]),
        )
        if (
            (row[0] < row[3])
            or (row[0] == row[3] and row[1] <= row[4])
            or (row[0] == row[3] and row[1] == row[4] and row[2] <= row[5])
        ):
            rows.append(row)
        else:
            rows.append(row[3:6] + row[0:3])

    sys.stderr.write(f"Removing duplicates from {len(rows)} pairs\n")
    rows = pd.DataFrame(rows)
    rows.sort_values(by=[0, 1, 3, 4], inplace=True)
    rows.drop_duplicates(inplace=True)

    sys.stderr.write(f"{rows.shape[0]} non-duplicate pairs to try merging\n")
    outfile.write("#chr1\tstart1\tend1\tchr2\tstart2\tend2\n")
    for (chr1, chr2), group in rows.groupby([0, 3]):
        sys.stderr.write("\r" + 80 * " ")
        sys.stderr.write(
            f"\rProcessing {group.shape[0]} pairs between {chrs[chr1]} and {chrs[chr2]}"
        )
        group = np.array(group)
        pairs = combos(np.array(group), args.dist)
        merged = new_rows(group, pairs)
        write_rows(merged, chrs, out=outfile, symetric=args.symetric)
    sys.stderr.write("\n")

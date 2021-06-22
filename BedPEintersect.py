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
def intersect(a, a0, a1, b, b0, b1):
    a, a0, a1 = row1
    b, b0, b1 = row2
    return (a == b) and ((a1 + dist) >= b0) and ((b1 + dist) >= a0)


@njit
def row_intersect(t1, t2, t3, q1, q2, q3, rtn):
    print(t1)
    for i in range(q1.shape):
        for j in range(t1.shape):
            if intersect(t1[j], t2[j], t3[j], q1[i], q2[i], q3[i]):
                rtn[i] = True
                break
    return rtn


def df_intersect(target, query):
    overlaps1 = row_intersect(
        target.chr.values,
        target.start.values,
        target.end.values,
        query.chr2.values,
        query.start2.values,
        query.end2.values,
        np.zeros(query.start2.values.shape, dtype=bool),
    )
    print(overlaps1)
    overlaps2 = row_intersect()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("target", help="positional input")
    parser.add_argument("query", help="positional input")
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
        "--cols",
        nargs="+",
        help="in order list of columns to use for second set of intervals",
        default=[3, 4, 5],
    )
    args = parser.parse_args()

    target = pd.read_csv(args.target, sep="\t").iloc[:, [0, 1, 2]]
    target.columns = ["chr", "start", "end"]
    query = pd.read_csv(args.query, sep="\t").iloc[
        :, [0, 1, 2, args.cols[0], args.cols[1], args.cols[2]]
    ]
    query.columns = ["chr", "start", "end", "chr2", "start2", "end2"]

    df_intersect(target, query)

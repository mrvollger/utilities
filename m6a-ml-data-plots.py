#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Docs here
"""

__author__ = "Mitchell R. Vollger"
__credits__ = ["Mitchell R. Vollger"]
__maintainer__ = "Mitchell R. Vollger"
__email__ = "mrvollger@gmail.com"
__status__ = "Development"


import logging
import argparse
import sys
import tqdm
import numpy as np
import logging
import pandas as pd
import os
import pickle
import gnuplotlib as gp


def plot(args):
    ml_data = np.load(args.input)
    windows = ml_data["features"]
    labels = ml_data["labels"]
    strands = ml_data["strands"]

    if args.n is not None and args.n < len(labels):
        logging.info(f"Subsampling to {args.n} data points")
        idx = np.random.choice(len(labels), args.n, replace=False)
        windows = windows[idx]
        labels = labels[idx]
        strands = strands[idx]

    logging.info(f"Data shape: {windows.shape}; Labels shape: {labels.shape}")
    logging.info(
        f"Positives: {labels.sum():,}\tNegatives: {len(labels) - labels.sum():,}"
    )

    for strand in [0, 1]:
        central_ip = (
            255 * windows[labels & (strand == strands), 4, args.window_size // 2]
        )
        non_m6a = 255 * windows[~labels & (strand == strands), 4, args.window_size // 2]
        p_base_ave = 255 * windows[labels & (strand == strands), 0:4, 5:10].mean(axis=0)
        n_base_ave = 255 * windows[~labels & (strand == strands), 0:4, 5:10].mean(
            axis=0
        )
        logging.info(f"Strand: {strand}")

        logging.info(
            f"Mean IPD at m6A:\t{central_ip.mean():.4g} +/- {central_ip.std():.4g}"
        )
        logging.info(f"Base weights:\n{p_base_ave}")
        logging.info(f"Base weights:\n{n_base_ave}")
        logging.info(
            f"Mean IPD at non-m6A:\t{non_m6a.mean():.4g} +/- {non_m6a.std():.4g}"
        )
        logging.info("")
        pos = np.mean(255 * windows[labels & (strand == strands), 4, :], axis=0)
        neg = np.mean(255 * windows[~labels & (strand == strands), 4, :], axis=0)
        gp.plot(
            pos,
            _with="lines",
            terminal="dumb 80,30",
            unset="grid",
            title=f"IPD at m6A on the {strand}",
        )
        gp.plot(
            neg,
            _with="lines",
            terminal="dumb 80,30",
            unset="grid",
            title=f"IPD at non-m6A on the {strand}",
        )


def parse():
    """Change all reads (and header) in a bam file to have one read group (RG)"""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input",
        help="Input npz file",
    )
    parser.add_argument(
        "-t", "--threads", help="Number of threads to use", type=int, default=8
    )
    parser.add_argument("-w", "--window-size", help="window size", type=int, default=15)
    parser.add_argument("-n", help="Max number of data points", type=int, default=None)
    parser.add_argument(
        "-v", "--verbose", help="increase logging verbosity", action="store_true"
    )
    args = parser.parse_args()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(format=log_format, level=log_level)
    return args


def main():
    args = parse()
    plot(args)


if __name__ == "__main__":
    main()

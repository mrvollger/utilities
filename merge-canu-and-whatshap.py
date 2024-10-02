#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Docs here"""

__author__ = "Mitchell R. Vollger"
__credits__ = ["Mitchell R. Vollger"]
__maintainer__ = "Mitchell R. Vollger"
__email__ = "mrvollger@gmail.com"
__status__ = "Development"

MAX_ROWS = 100

import logging
import argparse
import sys
import pandas as pd


def read_canu(file, tag):
    df = pd.read_csv(file, sep="\t", names=["seq"])
    df["canu_hap"] = tag
    df["seq"] = df.seq.str.strip(">")
    df.drop_duplicates(inplace=True)
    logging.info(f"Read {len(df):,} sequences from {file}")
    return df


def read_whatshap(file):
    df = pd.read_csv(
        file, sep="\t", names=["seq", "whatshap_hap", "phaseblock", "chr"], comment="#"
    )
    df.drop_duplicates(inplace=True)
    logging.info(f"Read {len(df):,} sequences from whatshap")
    logging.info(f"Read {len(df.phaseblock.unique()):,} phaseblocks from whatshap")
    logging.info(f"Whatshap hap counts:\n{df.whatshap_hap.value_counts()}")
    return df


def parse():
    """Change all reads (and header) in a bam file to have one read group (RG)"""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "canu",
        nargs=3,
        help="Input files, in order: canu paternal, canu maternal, canu unknown",
    )
    parser.add_argument(
        "whatshap",
        help="Input file",
    )
    parser.add_argument(
        "-c",
        "--prioritize-canu",
        help="keep phased canu reads even if they disagree with whatshap",
        action="store_true",
    )
    parser.add_argument(
        "-o", "--output-prefix", help="Output prefix", default="canu-plus-whatshap"
    )
    parser.add_argument(
        "-t", "--threads", help="Number of threads to use", type=int, default=8
    )
    parser.add_argument(
        "-v", "--verbose", help="increase logging verbosity", action="store_true"
    )
    args = parser.parse_args()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(format=log_format, level=log_level)
    return args


def assign_based_on_canu(group_df):
    canu_counts = group_df.canu_hap.value_counts()
    # logging.info(f"{canu_counts}")
    maternal = canu_counts.get("maternal", 0)
    paternal = canu_counts.get("paternal", 0)
    if paternal > maternal:
        group_df.merged_hap = "paternal"
    elif maternal > 0:
        group_df.merged_hap = "maternal"
    else:
        group_df.merged_hap = "unknown"
    return group_df


def main():
    args = parse()
    canu_dfs = []
    for canu, tag in zip(args.canu, ["paternal", "maternal", "unknown"]):
        canu_dfs.append(read_canu(canu, tag))
    canu_df = pd.concat(canu_dfs)
    whatshap_df = read_whatshap(args.whatshap)
    merged_df = whatshap_df.merge(canu_df, on="seq", how="left")

    merged_df["merged_hap"] = "unknown"
    merged_df = (
        merged_df[merged_df.whatshap_hap != "none"]
        .groupby(["phaseblock", "whatshap_hap"])
        .apply(assign_based_on_canu)
    )
    # count disagreements
    has_canu = merged_df.canu_hap != "unknown"
    disagreements = merged_df[has_canu].canu_hap != merged_df[has_canu].merged_hap
    logging.info(
        f"{disagreements.sum()/merged_df[has_canu].shape[0]:.2%} of whatshap calls disagree with the phased canu calls."
    )
    if args.prioritize_canu:
        logging.info("Prioritizing canu phasing.")

    # say # of reassignments
    phasing_added = (merged_df[~has_canu].merged_hap != "unknown").sum()
    logging.info(
        f"{phasing_added:,} unassigned sequences were assigned a haplotype based on whatshap phasing"
    )

    # assign final haplotype
    # logging.info(f"Canu counts:\n{merged_df.canu_hap.value_counts()}")
    merged_df["hap"] = merged_df.canu_hap
    merged_df.loc[~has_canu, "hap"] = merged_df.merged_hap[~has_canu]
    # logging.info(f"Merged counts:\n{merged_df.hap.value_counts()}")

    # make final outputs
    out = canu_df.merge(merged_df[["seq", "merged_hap", "hap"]], on="seq", how="left")
    out.loc[out.hap.isna(), "hap"] = out.canu_hap[out.hap.isna()]

    # drop ambiguous reads from phasing
    if not args.prioritize_canu:
        known = (out.canu_hap.isin(["maternal", "paternal"])) & (
            out.merged_hap.isin(["maternal", "paternal"])
        )
        disagreements = known & (out.canu_hap != out.merged_hap)
        out.loc[disagreements, "hap"] = "unknown"
        logging.info(
            f"{disagreements.sum()/out.shape[0]:.2%} of total reads changed to unknown due to whatshap and canu disagreements"
        )

    logging.info(f"Canu counts:\n{out.canu_hap.value_counts()}")
    logging.info(f"Final merged counts:\n{out.hap.value_counts()}")
    z = (out.canu_hap != "unknown").sum()
    logging.info(f"Canu phasing rate: {z/len(out):.2%}")
    z = (out.hap != "unknown").sum()
    logging.info(f"Merged phasing rate: {z/len(out):.2%}")

    for hap in ["paternal", "maternal", "unknown"]:
        out[out.hap == hap].seq.to_csv(
            f"{args.output_prefix}-{hap}.txt", index=False, header=False
        )


if __name__ == "__main__":
    main()

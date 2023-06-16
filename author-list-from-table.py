#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
import numpy as np

# global var for inputs
AFF = "Institution"
FIRST = "First Name"
MIDDLE = "Middle Name(s)/Initial(s)"
LAST = "Last Name"
EMAIL = "Email"
CORRESPONDING = "Corresponding Author"
FIRSTS = "Equal contribution"


def add_symbol(x, value="†"):
    if x == "" or x is None or pd.isna(x):
        return ""
    return value


def add_first(x):
    return add_symbol(x, value="*")


def add_corresponding(x):
    return add_symbol(x, value="†")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="positional input")
    parser.add_argument("-s", "--string", help="string option")
    parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
    parser.add_argument(
        "-l", "--list", nargs="*", help="list with zero or more entries"
    )
    parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    args = parser.parse_args()
    df = pd.read_csv(args.infile, sep="\t")
    df[MIDDLE] = df[MIDDLE].fillna("")

    df["AFFS"] = df[AFF].str.split(";")
    aff_nums = {}
    for aff in df[AFF]:
        affs = aff.split(";")
        for aff in affs:
            if aff not in aff_nums:
                aff_nums[aff] = len(aff_nums) + 1
    df["nums"] = df[AFF].apply(
        lambda x: ",".join([f"{aff_nums[aff]}" for aff in x.split(";")])
    )
    # first authors =

    df["Name"] = (
        (
            df[FIRST]
            + " "
            + df[MIDDLE]
            + " "
            + df[LAST]
            + df["nums"]
            + ","
            + df[FIRSTS].apply(add_first)
            + ","
            + df[CORRESPONDING].apply(add_corresponding)
        )
        .str.strip(",")
        .str.replace("  ", " ")
        .str.replace(",,", ",")
    )
    # print(df)

    print(", ".join(df.Name))
    print("")
    for aff, num in aff_nums.items():
        print(f"{num}. {aff.strip()}")
    print("\n† Corresponding author(s).")
    print("* These authors contributed equally to this work.")

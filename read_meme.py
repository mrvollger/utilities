#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import os
import sys
import argparse
from Bio.motifs import meme
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="positional input")
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    args = parser.parse_args()

    header = "Sequence name            Strand  Start   P-value                 Site    "
    breaker = "----------------------------------------------------------------"
    seen = False
    motif = ""
    with open(args.infile) as f:
        for line in f:
            if "Multilevel" in line:
                motif = line.strip().split()[1]
            elif breaker in line:
                seen = False
            elif header in line:
                seen = True
            elif seen and line[0] != "-":
                rec, strand, start, pval = line.strip().split()[0:4]
                offset = int(start)
                contig, st, en = re.match("(.+):([0-9]+)-([0-9]+)$", rec).groups()
                st, en = int(st), int(st) + offset
                print(f"{contig}\t{st}\t{en}\t{motif}")

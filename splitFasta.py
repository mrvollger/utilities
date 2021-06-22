#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("fasta", help=" input fasta")
parser.add_argument("out", help="output dir")
args = parser.parse_args()

import os
import sys
from Bio import SeqIO


recs = SeqIO.parse(args.fasta, "fasta")  # load the input fasta

for rec in recs:  # iterate through multi fasta
    newfile = os.path.abspath(
        args.out + "/" + rec.id + ".fasta"
    )  # make the name of the output file
    SeqIO.write(
        rec, newfile, "fasta"
    )  # write this record "rec" to an output file "newfile" in this "fasta" format

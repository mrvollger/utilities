#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument(
    "fasta",
    nargs="?",
    help="fasta file to mask",
    type=argparse.FileType("r"),
    default=sys.stdin,
)
parser.add_argument(
    "bed",
    nargs="?",
    help="bed file with regions to mask",
    type=argparse.FileType("r"),
    default=sys.stdin,
)
parser.add_argument(
    "-o",
    "--out",
    help="output fasta file",
    type=argparse.FileType("w"),
    default=sys.stdout,
)
parser.add_argument(
    "-a",
    "--add",
    help="a second fasta file to add to the output",
    type=argparse.FileType("r"),
    default=None,
)
parser.add_argument("-d", action="store_true", default=False)
args = parser.parse_args()

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO


def bed_to_dict(bed):
    mask = {}
    for line in bed:
        contig, start, end = line.strip().split()[0:3]
        start = int(start)
        end = int(end)
        if contig not in mask:
            mask[contig] = []
        mask[contig].append((start, end))
    return mask


mask = bed_to_dict(args.bed)
recs = SeqIO.parse(args.fasta, "fasta")

for rec in recs:
    sys.stderr.write("{}\n".format(rec.name))
    if rec.id in mask:
        seq = str(rec.seq)
        for start, end in mask[rec.id]:
            seq = seq[0:start] + "N" * (end - start) + seq[end:]
        rec.seq = Seq(seq, generic_dna)

    SeqIO.write(rec, args.out, "fasta")

if args.add is not None:
    recs = SeqIO.parse(args.add, "fasta")
    SeqIO.write(recs, args.out, "fasta")

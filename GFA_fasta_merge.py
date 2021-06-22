#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("GFA", nargs="?", help="input GFA")
parser.add_argument("fasta", nargs="?", help="input fasta")
parser.add_argument(
    "outfile",
    nargs="?",
    help="output bam file",
    type=argparse.FileType("w"),
    default=sys.stdout,
)
parser.add_argument("-d", action="store_true", default=False)
args = parser.parse_args()

from Bio import SeqIO

recs_d = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))


for line in open(args.GFA):
    tokens = line.strip().split("\t")
    if tokens[0] == "S" and tokens[2] == "*":
        rec_id = tokens[1]
        seq = str(recs_d[rec_id].seq)
        end = "\t".join(tokens[3:])
        args.outfile.write("{}\t{}\t{}\t{}\n".format("S", rec_id, seq, end))
    else:
        args.outfile.write(line)

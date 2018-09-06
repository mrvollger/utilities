#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("fasta", help="input fasta file" )
parser.add_argument("name",help="fasta name to grab, or file to list of names")
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()

from Bio import SeqIO
import os

toFind = []
if(os.path.exists(args.name)):
	for line in open(args.name):
		toFind.append(line.strip())
else:
	toFind.append(args.name)








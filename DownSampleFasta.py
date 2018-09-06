#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument('-i', '--infasta', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o','--outfasta', type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument("--perCutoff", 
		help="percentile cutoff, an input of 90 would return the top 10 percentile of reads", default=None )
parser.add_argument("--kbCutoff", type = float, default = None, help="kb cutoff only use reads of length [input] kb and greater")
args = parser.parse_args()

from Bio import SeqIO
import numpy as np
Cutoff = 1000000000000000000000000000000
# input setup
reads_dict = SeqIO.to_dict(SeqIO.parse(args.infasta, "fasta"))
length_dict = {}

def getLengths():
	global length_dict
	for readid in reads_dict:
		read = reads_dict[readid]
		length_dict[readid] = len(read.seq)

def keysByPerCutoff():
	global Cutoff
	if(args.perCutoff is None):
		return(set(reads_dict.keys()))
	# determine percential cutoff 	
	npa = np.array( list(length_dict.values()) )
	cutoff = np.percentile(npa, args.perCutoff)
	Cutoff = cutoff
	# determine reads to keep	
	tokeep = []
	for readid in length_dict:
		length = length_dict[readid]
		if(length >= cutoff):
			tokeep.append(readid)
	return(set(tokeep))

def keysByKbCutoff():
	global Cutoff 
	if(args.kbCutoff is None):
		return(set(reads_dict.keys()))
	kbCut = args.kbCutoff * 1000
	if(kbCut > Cutoff):
		Cutoff = kbCut
	tokeep = []
	for readid in length_dict:
		length = length_dict[readid]
		if(length >= kbCut):
			tokeep.append(readid)
	return( set(tokeep))


getLengths()
perkeys = keysByPerCutoff()
kbkeys = keysByKbCutoff()

print("THE USED CUTOFF WAS: {}".format( int(Cutoff )), file = sys.stderr)

towrite = perkeys.intersection(kbkeys)

for key in sorted(towrite):
	SeqIO.write(reads_dict[key], args.outfasta, "fasta")




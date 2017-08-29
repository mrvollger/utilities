#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="extract bam entries based on a list of read names")
parser.add_argument("bam", help="the bam file" )
parser.add_argument("names",help="list of read names",default=None)
parser.add_argument("--out", help="the name of out output bam file", default="extracted.bam" )
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
DEBUG=args.d

import pysam 


names = open(args.names).read().splitlines()
numNames=len(names)

bam = pysam.AlignmentFile(args.bam)
bamout= pysam.AlignmentFile(args.out, 'wb', header=bam.header)

counter = 0 
for record in bam.fetch(until_eof=True):
    if(record.query_name in names):
        counter += 1
        bamout.write(record)

#print(counter, numNames)
assert counter == numNames


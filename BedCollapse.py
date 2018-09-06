#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="This is similar to merge in bedtools, but instead of mergeing overlapping regions, it collapses all regions that are entirly contained within another region")
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()


import intervaltree as it
import pandas as pd
from sets import Set


bed = pd.read_csv( args.infile, header=None, delimiter=r"\s+")
numCols = len(list(bed))

# add all things to a data column 
bed['data'] = bed[  range(3,numCols)  ].apply(lambda x: '\t'.join(x.map(str))[1:], axis=1)
bed['key'] = bed[0] + "_" + bed[1].map(str) + "_" + bed[2].map(str)


keepKeys = Set([])
chrs = bed.groupby([0])
Safe = False # 99% sure it is safe to set Safe to false, and it is much faster
for CHR, regions in chrs:
	tuples = zip(regions[1], regions[2])
	tree  = it.IntervalTree.from_tuples(tuples)	
	contained = it.IntervalTree()

	for interval in tree.items():
		start = interval.begin 
		end = interval.end
		data = interval.data
		#tree.remove_envelop(start, end) 
		
		toRemove = tree.search(start, end, strict=True)
		size1 = len(toRemove)
			#if( size1 != 0):
		toRemove.discard( interval )
			#size2 = len(toRemove)
			#assert size1 == size2 + 1, "{}  {}".format(size1, size2)
		
		if(Safe): # I know this way is safe because I am not removing from my iterable tree
			contained.update(toRemove)	
		else: # might not be safe becuase I am removing form my iterator 	
			for rmMe in toRemove:
				tree.discard(rmMe)

		# the remove command deletes itself so I need to add it back in

	if(Safe):
		final = tree - contained
	else:
		final = tree
	
	for inv in final.items():
		keepKeys.add( "{}_{}_{}".format(CHR, inv.begin, inv.end))
	
#
bedOut = bed.loc[ bed['key'].isin(keepKeys) ]
# sort like a bed file and then on the last column in the bed file
bedOut.sort_values( [0,1, numCols-1], inplace=True)
bedOut.drop_duplicates(subset = [0,1,2], keep="first", inplace=True)

outCols = list(  Set(bed).difference(Set(["data", "key"])) )
bedOut[ outCols ].to_csv(args.outfile, sep="\t", header = False, index=False)



#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("sam",help="output png, the name must have the .png ext",default=None)
parser.add_argument("png",help="output png, the name must have the .png ext",default=None)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
DEBUG=args.d

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.rc('font', family='serif')
import pysam 
import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


samfile = pysam.AlignmentFile(args.sam)
contigs = []
lengths = {}
for read in samfile.fetch():
	if(read.reference_name not in contigs):
		lengths[read.reference_name] = read.reference_length
		contigs.append(read.reference_name)
samfile.close()

numGroups = len(contigs)
fig, axsList = plt.subplots(numGroups, figsize=(16,2*numGroups) )
# convert to array
if(type(axsList) != type(np.array([0]))):
    axsList = np.array([axsList])

# convert to dictionary 
axs={}
for idx, ax in enumerate(axsList):
	# associate a ax with a specific contig
	axs[ contigs[idx] ] = ax
	axs[ contigs[idx] ].set_xlim([0, lengths[contigs[idx]]])
	axs[ contigs[idx] ].set_ylim([0, 1000])


def addPatchToAx(ax, x, y, width, height, operation):
	col = {2:"gray", 7:"darkgreen", 8:"darkred"}
	ax.add_patch( patches.Rectangle((x,y), width, height, alpha=0.75, linewidth=0, color = col[operation] )) 


def squaresFromCigar(ax, cigar, start, y, height):
	x = start 
	for pair in cigar:
		operation = pair[0]
		width = pair[1]
		if(operation in [2,7,8]):
			addPatchToAx(ax, x, y, width, height, operation)
			x += width
	return(x)
def getPerID(aln):
	stats = aln.get_cigar_stats()[0]
	match = stats[7]*1.0
	mismatch = stats[8]*1.0
	ins = stats[1]*1.0
	dele = stats[2]*1.0
	rtn = "{0:.2f}".format( match/(match + mismatch + ins + dele)*100  ) 
	print(rtn)
	return(rtn)

# add contigs from sam file
if(args.sam is not None):
	samfile = pysam.AlignmentFile(args.sam)
	# a map from reference to names of queriers that mathc  
	ref={}
	reads={}
	for read in samfile.fetch():
		if(read.is_unmapped):
			continue 
		# add reads to a dictionary by query
		reads[read.query_name] = read
		# add queries to a dictionary by contig they map to
		if(read.reference_name not in ref):
			ref[read.reference_name] = []
		if(read.reference_name != read.query_name):
			ref[read.reference_name].append(read.query_name)

	for contig in sorted(ref):
		print("ploting " + str(contig) )
		ax = axs[contig]
		queries = ref[contig]
		if(len(queries) == 0):
			continue 
		max_y = 1000 #max(group["depth"])
		height = max_y/(len(queries))
		# add rectangles for each alignment / cigar string tupple 
		for counter, query in enumerate(queries):
			y = counter * max_y/(len(queries))
			read = reads[query]
			end = squaresFromCigar(ax, read.cigar, read.reference_start, y , height)
			assert(end == read.reference_end)
			perID= getPerID(read)
			label="{}:{}-{};{};{};{}".format(query, read.query_alignment_start, read.query_alignment_end, 
					read.is_reverse, read.infer_query_length(), perID)
			ax.text(read.reference_start, y, label)
		ax.set_title(contig)


# this shoudl fix overlapping lables 
plt.tight_layout()
plt.show()
plt.savefig(args.png)


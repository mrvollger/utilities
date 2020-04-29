#!/usr/bin/env python
import argparse
import os 
import sys
import numpy as np

# global var for inputs
args=None 

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("infile", help="positional input")
	parser.add_argument("-c", "--column", help="column to read from", type=int, default=1)
	parser.add_argument("-b", "--bins", help="number of bins to make", type=int, default=30)
	parser.add_argument("-w", "--width", help="number of bins to make", type=int, default=80)
	parser.add_argument('-l', "--log", help="use log10 base bins",  action="store_true", default=False)
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()

	data = []
	for line in open(args.infile):
		data.append( float(line.strip().split()[args.column -1 ]) )
	data = np.array(data)
	MAX = data.max()
	MIN = data.min()
	bins = np.linspace(MIN, MAX, args.bins)
	if(args.log):
		LMIN = np.floor(np.log10(MIN))
		LMAX = np.ceil(np.log10(MAX))
		bins=np.logspace(LMIN, LMAX, LMAX-LMIN+1)


	inds = np.digitize(data, bins)
	maxbincount = np.max(np.bincount(inds))
	
	#print(inds, inds.shape, data.shape, bins.shape, inds.min(), inds.max(), maxbincount)
	sys.stdout.write("[start\tend)\tcount\thistogram\n")
	
	for b in range( bins.shape[0] - 1):
		count = (inds == b+1).sum()
		n = count / maxbincount * args.width 
		text = "*"*int(n)
		if(int(n) -n != 0): text += "."
		sys.stdout.write(f"[{bins[b]:.2f}\t{bins[b+1]:.2f})\t{count}\t{text}\n")




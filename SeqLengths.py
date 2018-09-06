#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="Get the lengths of reads from fasta/fastq/(.gz) files.")
parser.add_argument("infiles", nargs="*", help="input fasta/fastq file")
parser.add_argument("-f", "--fofn", help="input FOFN of fasta/fastq files, replaces infiles argument")
parser.add_argument("-o", "--out", help="output lengths file. Has the form: readname, length", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-p',"--plot", default=None,help="Specify a place to plot a hsitogram of read lengths")
parser.add_argument("--plotwith", default=None,help="If you have already calcualted readlength before you can specify plotwith to avoid rereading the fasta files")
parser.add_argument('-t', "--threads", type=int, default=1, help="number of threads to use, only used in reading in fasta/fastq files")
args = parser.parse_args()


from Bio import SeqIO
import gzip
import mimetypes
from multiprocessing.dummy import Pool as ThreadPool 
pool = ThreadPool(args.threads)



def readSeqs(f, ftype):
	lengths = []
	for rec in SeqIO.parse(f, ftype):
		lengths.append( (rec.name, len(rec.seq)) )
	return(lengths)

#for seqfile in args.infiles:
def readFile(seqfile):
	#mime = mimetypes.guess_type(seqfile)
	#print(mime)
	if(seqfile[-3:] == ".gz"):
		f = gzip.open(seqfile, 'rt')
	else:
		f = open(seqfile)
	
	# get file type 
	ftype = None
	tags = seqfile.split(".")
	if("fasta" in tags):
		ftype = "fasta"
	elif("fastq" in tags):
		ftype = "fastq"
	elif("fq" in tags):
		ftype = "fastq"
	elif("fa" in tags):
		ftype = "fasta"
	assert ftype is not None, "File type unknown"
	print("File:{}\tFile Type:{}".format(seqfile, ftype))

	lengths = readSeqs(f, ftype)
	f.close()
	return(lengths)

def NX(nums, X):
	nums = sorted(nums, key=int, reverse=True)
	datathresh = sum(nums)*(X/100.0)	
	total = 0
	for num in nums:
		total += num
		if(total >= datathresh):
			return(num)
	return(0)



flat_list=None
if(args.plotwith is None):
	myfiles = args.infiles
	if(args.fofn is not None):
		myfiles = []
		for line in open(args.fofn).readlines():
			myfiles.append(line.strip())
	lengths = pool.map(readFile, myfiles)

	flat_list = []
	for l in lengths:
		flat_list += l

	out = ""
	for name, length in flat_list:
		out += "{}\t{}\n".format(name, length)
	args.out.write(out)
	plotwith = args.out


if(args.plot is not None):
	import pandas as pd
	import numpy as np
	import matplotlib
	matplotlib.use('Agg')
	import seaborn as sns
	import matplotlib.pyplot as plt
	sns.set(font_scale=2)
	sns.set_style("ticks")
	fig, ax = plt.subplots(figsize=(16,9))

	if(args.plotwith is not None):
		df = pd.read_csv(plotwith, sep="\t", header=None, names=["name", "length"])
	elif(flat_list is not None):
		df = pd.DataFrame(flat_list, columns=["name", "length"])
	else:
		print("No data to plot with, try --plotwith", file = sys.stderr)

	df["LogLength"]=np.log10( np.clip(df["length"],10,None) )

	# make histogram	
	sns.distplot(df.LogLength, bins=100, kde=False, rug=False, ax=ax)

	# get maxes
	myy = plt.gca().get_ylim()[1]
	mymax = max(df["length"])
	
	# add vertical lines
	vals = [df["length"].median(), df["length"].mean(), 
			NX(df["length"], 50.0), NX(df["length"], 1.0), mymax]
	names = ["Median", "Mean", "N50", "N1", "Max"]
	divs = [.9, .8, .7, .6, .5]
	for name, val, div in zip(names, vals, divs):
		plt.axvline(x=np.log10(val), color="darkred", linestyle="dashed")
		label = "{}={}".format(name, round(val/1000, 1))
		plt.text(np.log10(val)+.01, myy*div, label, fontsize=16)
	
	# make x ticks 	
	xts = [2,3,4,5,6,7]
	xls = ["0.1","1", "10", "100", "1,000", "10,000"]
	if(np.log10(mymax) <= 6):
		xts = xts[:-1]; xls=xls[:-1]

	# plot 
	plt.xticks(xts, xls)
	Gb = sum(df["length"])/1000000000
	plt.xlabel("Length (kbp), Total Gb {:.2f}".format(Gb))
	plt.savefig(args.plot)








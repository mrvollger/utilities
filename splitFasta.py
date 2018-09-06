#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("fasta", help=" input fasta" )
parser.add_argument("out",help="output dir",default="reference_fractionated")
args = parser.parse_args()

import glob
import os
import sys
import re
from Bio import SeqIO
import pysam 
import runCmd
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


print("Staritng")
recs = list(SeqIO.parse(args.fasta, "fasta"))
numdig=len(str(len(recs)))
form = args.out  + '/{}_{' + ":0{}".format(numdig) + '}'
print(form)
conversion = ""
for idx, rec in enumerate(recs):
	myfile = form.format(os.path.basename(args.fasta), idx)
	newname = "contig{}".format(idx)
	conversion += "{}\t{}\n".format(newname, rec.id)
	rec.id = newname 
	rec.name = newname
	rec.description = ""
	SeqIO.write(rec, myfile, "fasta")


open("namechange.txt", "w+").write(conversion)

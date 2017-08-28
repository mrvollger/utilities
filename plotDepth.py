#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("depth", help="tab seperated file, contig   pos   depth" )
parser.add_argument("png",help="output png, the name must have the .png ext",default=None)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
DEBUG=args.d

import glob
import os
import sys
import re
import itertools 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam 
import runCmd
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


colNames = ["contig", "pos", "depth"]
df = pd.read_csv(args.depth, sep="\t", names=colNames)

#df.plot(x="pos", y="depth")
grouped = df.groupby('contig')
fig, axs = plt.subplots(len(grouped))
counter = 0
for key, group in grouped:
    print(key)
    group.plot(ax=axs[counter], kind='scatter', x='pos', y='depth', label=key)
    counter += 1
plt.show()
plt.savefig(args.png)


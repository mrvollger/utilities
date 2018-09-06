#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input bam file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outfile",nargs="?", help="output bam file", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()

import glob
import os
import sys
import re
import itertools 
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam 
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO





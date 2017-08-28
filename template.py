#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("", help="" )
parser.add_argument("",help="",default=None)
parser.add_argument("--", help="", default=None )
parser.add_argument("--", help="", default="" )
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
 = args.
 = args.
 = args.
 = args.
DEBUG=args.d

import glob
import os
import sys
import re
import itertools 
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam 
import runCmd
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO





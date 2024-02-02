#!/usr/bin/env python
import argparse
import multiprocessing
import sys
from multiprocessing import current_process

import pysam

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input bam file")
parser.add_argument("--prefix", default="chunk")
parser.add_argument("-t", "--threads", help="threads to use", type=int, default=1)
parser.add_argument("-d", action="store_true", default=False)
args = parser.parse_args()

infile = pysam.AlignmentFile(args.infile)
outfs = {
    "ForkPoolWorker-{}".format(idx): pysam.AlignmentFile(
        "{}.{}.bam".format(args.prefix, idx), "wb", template=infile
    )
    for idx in range(1, args.threads + 1)
}


def write_rec(rec):
    name = current_process().name
    outf = outfs[name]
    outf.write(rec)
    return rec


def bam_writer():
    counter = 0
    pool = multiprocessing.Pool(args.threads)
    bam_iter = infile.fetch(until_eof=True)
    for rec in pool.imap_unordered(write_rec, bam_iter, chunksize=100):
        counter += 1
        if counter % 100 == 0:
            sys.stderr.write("\rWritten {}".format(counter))


print("here", __name__)
if __name__ == "__main__":
    print("hello")
    bam_writer()

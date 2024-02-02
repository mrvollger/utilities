#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("bam", help="input files (bam, sam, cram)")
parser.add_argument(
    "-r", "--regions", nargs="*", help="region(s) in contig:start-end format"
)
parser.add_argument(
    "-o",
    "--outfile",
    help="output file (sam), default is stdout",
    type=argparse.FileType("w"),
    default=sys.stdout,
)
parser.add_argument("--bed", help="bed file with regions", default=None)
args = parser.parse_args()

import os
import re
import pysam
import collections
import subprocess
import tempfile

M = 0  # M  BAM_CMATCH      0
I = 1  # I  BAM_CINS        1
D = 2  # D  BAM_CDEL        2
N = 3  # N  BAM_CREF_SKIP   3
S = 4  # S  BAM_CSOFT_CLIP  4
H = 5  # H  BAM_CHARD_CLIP  5
P = 6  # P  BAM_CPAD        6
E = 7  # =  BAM_CEQUAL      7
X = 8  # X  BAM_CDIFF       8
B = 9  # B  BAM_CBACK       9
NM = 10  # NM       NM tag  10
conRef = [M, D, N, E, E]  # these ones "consume" the reference
conQuery = [M, I, S, E, X]  # these ones "consume" the query
conAln = [M, I, D, N, S, E, X]  # these ones "consume" the alignments


def exists(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return True
    else:
        return False


def get_regions():
    regions = []
    if args.regions is not None:
        for region in args.regions:
            match = re.match("(.*):(\d*)-(\d*)", region)
            if match is None:
                print("WARNING: no match for " + region, file=sys.stderr)
            else:
                regions.append(
                    (match.group(1), int(match.group(2)) - 1, int(match.group(3)))
                )

    if args.bed is not None:
        for region in open(args.bed):
            region = region.strip().split()
            regions.append((region[0], int(region[1]), int(region[2])))

    assert len(regions) > 0, "ERROR: no regions specified!"
    return regions


def rle(opts):
    # print(opts, len(opts))
    pre = None
    # add a terminating opt so it returns the last start of rle
    opts += [None]
    count = 1
    rtn = []
    idx = 0
    for cur in opts:
        if cur != pre and pre is not None:
            # if cur is none we are at tehn end and need to adjsut our counts by 1
            if cur is None:
                count -= 1

            rtn.append((pre, count))
            count = 1
        else:
            count += 1
        pre = cur

    return rtn


def change_alnseg(read, all_opts, q_st, q_en, r_st, r_en):
    # subset sequence
    read.seq = read.seq[q_st:q_en]

    # susbset quality values
    if read.query_alignment_qualities:
        read.query_alignment_qualities = read.query_alignment_qualities[q_st:q_en]

    # change alignment start location
    read.reference_start = r_st

    # remove tags
    read.set_tags(None)

    # change cigar
    read.cigartuples = rle(all_opts)
    # print(rle(all_opts))

    return read


def trim_read(read, contig, start, end):
    q_st = 0
    r_st = 0
    q_en = 0
    r_en = 0

    rpos = read.reference_start
    qpos = read.query_alignment_start
    all_opts = []
    in_rgn = False
    for opt, length in read.cigartuples:
        for i in range(length):
            # check for start adn end poses
            if rpos == start:
                q_st = qpos
                r_st = rpos
                in_rgn = True
            if rpos == end:
                q_en = qpos
                r_en = rpos
                in_rgn = False

            # incremend position
            if opt in conRef:
                rpos += 1
            if opt in conQuery:
                qpos += 1

            if in_rgn:
                all_opts.append(opt)
            else:
                all_opts.append(H)

    return change_alnseg(read, all_opts, q_st, q_en, r_st, r_en)


if __name__ == "__main__":
    inbam = pysam.AlignmentFile(args.bam)
    out = pysam.AlignmentFile(args.outfile, "w", template=inbam)
    counter = 0
    regions = get_regions()
    for contig, start, end in regions:
        for read in inbam.fetch(contig=contig, start=start, stop=end):
            read = trim_read(read, contig, start, end)
            sys.stderr.write("\rReads Proccesed: {}".format(counter))
            counter += 1
            out.write(read)
            break

    sys.stderr.write("\n")
    out.close()

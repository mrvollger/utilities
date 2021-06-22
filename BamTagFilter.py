#!/usr/bin/env python
import argparse
import sys
import os

parser = argparse.ArgumentParser(
    description="extract bam entries based on a list of read names"
)
parser.add_argument(
    "-t",
    "--tags",
    nargs="+",
    help='space seperated list of tags opteration value in quotes, e.g. "NM>10" "RG=FAKE"',
)
parser.add_argument("-b", "--bams", help="the bam file(s)", nargs="+")
parser.add_argument(
    "-o", "--out", help="out bam file", type=argparse.FileType("w"), default=sys.stdout
)
args = parser.parse_args()

import pysam
import re
import tempfile

# import collections
import subprocess


def merged_header(infiles):
    tmps = []
    notyet = []
    for bam in infiles:
        tmp = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        # print(tmp.name, file = sys.stderr)
        subprocess.call(["samtools", "view", "-H", "-b", bam], stdout=tmp)
        # for line in open(tmp.name):
        #       print(line)
        tmps.append(tmp.name)
        notyet.append(tmp)

    merged = tempfile.NamedTemporaryFile(suffix=".bam", delete=False).name
    # merged = "tmp.header.bam"
    subprocess.call(["samtools", "merge", "-f", merged] + tmps)
    bam = pysam.AlignmentFile(merged, check_sq=False)
    return bam


def hasTag(read, tag):
    # print(tag, file=sys.stderr)
    match = re.match(r"([A-Za-z][A-Za-z0-9])(==|>=|<=|>|<|=)(.+)", tag)
    assert match is not None, "Tag format is bad" + tag
    # print(match.groups(), file = sys.stderr)
    tag, opt, val = match.groups()

    result = read.get_tag(tag)
    mytype = type(result)
    # convert val to correct type
    val = mytype(val)
    if opt == ">":
        return result > val
    elif opt == "<":
        return result < val
    elif opt == "=":
        return result == val
    else:
        print("unpsorted opt: " + opt, file=sys.stderr)
        exit(1)


print(args.tags, file=sys.stderr)


super_header = merged_header(args.bams)
bamout = pysam.AlignmentFile(args.out, "wb", template=super_header)


counter = 0
counter2 = 0
for bamf in args.bams:
    print(bamf, file=sys.stderr)
    bam = pysam.AlignmentFile(bamf, check_sq=False)
    for read in bam.fetch(until_eof=True):
        keep = True
        for tag in args.tags:
            has = hasTag(read, tag)
            if has is False:
                keep = False
            # print(has, tag, file=sys.stderr)

        if keep:
            bamout.write(read)
            counter2 += 1
        counter += 1
        if counter % 100000 == 0:
            print(counter, counter2, file=sys.stderr)

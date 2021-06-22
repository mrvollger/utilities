#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument(
    "-r", "--regions", nargs="*", help="region(s) in contig:start-end format"
)
parser.add_argument("-i", "--infiles", nargs="+", help="input files (bam, sam, cram)")
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


def exists(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return True
    else:
        return False


def old_merged_header(infiles):
    super_header = collections.OrderedDict()
    bams = []
    for bamn in infiles:
        if exists(bamn):
            bam = pysam.AlignmentFile(bamn)
            header = bam.header.as_dict()

            for key in header:
                val = header[key]
                # append to super header
                if type(val) == type(dict()):
                    # add appropriate class to the super header
                    if key not in super_header:
                        super_header[key] = {}

                    for subkey in val:
                        subval = val[subkey]
                        if subkey in super_header[key]:
                            if subkey in ["SO", "GO"]:
                                assert super_header[key][subkey] == subval, (
                                    key,
                                    val,
                                    super_header[key],
                                )
                            if subkey == "VN" and super_header[key][subkey] != subval:
                                print(
                                    "WARNING: version conflict",
                                    super_header[key][subkey],
                                    subval,
                                    file=sys.stderr,
                                )
                        else:
                            super_header[key][subkey] = subval

                elif type(val) == type(list()):
                    # add appropriate class to the super header
                    if key not in super_header:
                        super_header[key] = []

                    for subval in val:
                        if subval not in super_header[key]:
                            super_header[key].append(subval)
            bams.append(bam)

    return (super_header, bams)


def merged_header(infiles):
    tmps = []
    notyet = []
    for bam in infiles:
        tmp = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        # print(tmp.name, file = sys.stderr)
        subprocess.call(["samtools", "view", "-H", "-b", bam], stdout=tmp)
        # for line in open(tmp.name):
        # 	print(line)
        tmps.append(tmp.name)
        notyet.append(tmp)

    merged = tempfile.NamedTemporaryFile(suffix=".bam", delete=False).name
    # merged = "tmp.header.bam"
    subprocess.call(["samtools", "merge", "-f", merged] + tmps)
    bam = pysam.AlignmentFile(merged, check_sq=False)
    return bam


def get_regions():
    regions = []
    if args.regions is not None:
        for region in args.regions:
            match = re.match("(.*):(\d*)-(\d*)", region)
            if match is None:
                print("WARNING: no match for " + region, file=sys.stderr)
            else:
                regions.append(
                    (match.group(1), int(match.group(2)), int(match.group(3)) + 1)
                )

    if args.bed is not None:
        for region in open(args.bed):
            region = region.strip().split()
            regions.append((region[0], int(region[1]), int(region[2])))

    assert len(regions) > 0, "ERROR: no regions specified!"
    print(regions, file=sys.stderr)
    return regions


regions = get_regions()

super_header = merged_header(args.infiles)

out = pysam.AlignmentFile(args.outfile, "w", template=super_header)
print("Fetching regions", file=sys.stderr)
for bamf in args.infiles:
    if exists(bamf):
        bam = pysam.AlignmentFile(bamf)
        for contig, start, end in regions:
            for read in bam.fetch(contig=contig, start=start, stop=end):
                out.write(read)
out.close()

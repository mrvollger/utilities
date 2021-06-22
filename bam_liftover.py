#!/usr/bin/env python
import argparse
import os
import sys
import pysam

# global var for inputs
args = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bed", help="positions to lift over")
    parser.add_argument("bam", help="bam/cram to process, - for stdin")
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)

    bed = [line.strip().split() for line in open(args.bed)]
    for rgn in bed:
        rgn[1] = int(rgn[1])
        rgn[2] = int(rgn[2])

    for rgn in bed:
        ranges = {}
        start = rgn[1]
        end = rgn[2]
        contig = rgn[0]
        for rec in bam.fetch(contig=contig, start=start, end=end):
            name = rec.query_name
            length = rec.infer_read_length()
            offset = 0
            if rec.cigartuples[0][0] == 5:  # is hard clipped
                offset = rec.cigartuples[0][1]
            rec.cigartuples
            qstart = None
            rstart = None
            qend = None
            rend = None
            for qp, rp in rec.get_aligned_pairs(matches_only=True):
                if rp < start or rp >= end:
                    continue
                if qstart is None or qp < qstart:
                    qstart = qp
                    rstart = rp
                elif qend is None or qp > qend:
                    qend = qp
                    rend = rp

            qstart += offset
            qend += offset
            if rec.is_reverse:
                tmp = qstart
                qstart = length - qend
                qend = length - tmp
            print(name, qstart, qend, contig, rstart, rend)

    print(bed)

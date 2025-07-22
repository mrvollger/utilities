#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional
import pysam


def main(
    infile: Optional[Path] = None,
    outfile: Optional[Path] = None,
    *,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    :param count: Number of times to display the greeting
    :param verbose: Set the logging level of the function
    """
    if infile is None:
        infile = sys.stdin
    if outfile is None:
        outfile = sys.stdout

    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    bam = pysam.AlignmentFile(infile)
    for rec in bam.fetch(until_eof=True):
        bp, count = rec.get_cigar_stats()
        i_bp = bp[1]
        i_count = count[1]
        d_bp = bp[2]
        d_count = count[2]

        print(f"{rec.query_name}\t{i_bp}\t{i_count}\t{d_bp}\t{d_count}")

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

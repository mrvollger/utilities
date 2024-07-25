#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional
import pysam
from tqdm import tqdm


def run(bam, nrows):
    A = ("A", 0, "a")
    T = ("T", 1, "a")
    counts = {A: 0, T: 0}
    read_counts = {A: 0, T: 0}
    for idx, rec in tqdm(enumerate(bam.fetch(until_eof=True))):
        if idx == nrows:
            break
        bm = rec.modified_bases_forward
        for key in counts:
            if key in bm:
                counts[key] += len(bm[key])
                read_counts[key] += 1

    print("Strand\tCount\tReadCount")
    print(f"A+a\t{counts[A]:,}\t{read_counts[A]:,}")
    print(f"T-a\t{counts[T]:,}\t{read_counts[T]:,}")


def main(
    infile: Optional[Path] = None,
    *,
    verbose: int = 0,
    threads: int = 16,
    nrows: int = 100_000,
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

    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    bam = pysam.AlignmentFile(infile, "rb", threads=threads)
    run(bam, nrows)
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

#!/usr/bin/env python
import logging
import sys
from pathlib import Path
from typing import Optional

import defopt
import pysam
from tqdm import tqdm


def run(bam, min_mapq=0):
    h1_count = 0
    h2_count = 0
    other_count = 0
    for rec in tqdm(bam.fetch(until_eof=True)):
        if rec.is_secondary or rec.is_supplementary:
            continue

        if rec.has_tag("HP"):
            hp = rec.get_tag("HP")
            if hp == 1:
                h1_count += 1
            elif hp == 2:
                h2_count += 1
        else:
            other_count += 1

    total = h1_count + h2_count + other_count
    phased = h1_count + h2_count
    print(f"{h1_count}, {h2_count}, {other_count}, {total}")
    print(f"{phased / total:%}")


def main(
    infile: Optional[Path] = None,
    *,
    min_mapq: int = 0,
    threads: int = 8,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    :param count: Number of times to display the greeting
    :param verbose: Set the logging level of the function
    """

    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    bam = pysam.AlignmentFile(infile, "rb", threads=threads)
    run(bam, min_mapq=min_mapq)

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

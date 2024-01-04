#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional
import pysam
from tqdm import tqdm


def run(bam, obam, hap1_tag="hap1", hap2_tag="hap2", min_mapq=0):
    for rec in tqdm(bam.fetch(until_eof=True)):
        if rec.is_unmapped or rec.mapping_quality <= min_mapq:
            rec.set_tag("HP", None)
        elif hap1_tag in rec.reference_name:
            rec.set_tag("HP", 1)
        elif hap2_tag in rec.reference_name:
            rec.set_tag("HP", 2)
        else:
            rec.set_tag("HP", None)
        obam.write(rec)


def main(
    infile: Optional[Path] = None,
    outfile: Optional[Path] = None,
    *,
    hap1_tag: str = "haplotype1",
    hap2_tag: str = "haplotype2",
    min_mapq: int = 1,
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
    if infile is None:
        infile = sys.stdin
    if outfile is None:
        outfile = sys.stdout

    bam = pysam.AlignmentFile(infile, "rb", threads=threads)
    obam = pysam.AlignmentFile(outfile, "wb", template=bam, threads=threads)
    run(bam, obam, min_mapq=min_mapq, hap1_tag=hap1_tag, hap2_tag=hap2_tag)

    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

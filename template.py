#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional


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

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

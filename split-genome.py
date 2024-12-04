#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path


def main(
    fai: Path,
    outfiles: list[Path],
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
    logging.basicConfig(level=1)
    n = len(outfiles)
    total_len = 0
    chroms = []
    with open(fai, "r") as f:
        for line in f:
            chrom, length, *_ = line.strip().split("\t")
            chroms.append((chrom, int(length)))
            total_len += int(length)

    outfiles = [open(f, "w") for f in outfiles]
    split_size = total_len // n + 1
    
    bases_will_be_covered = split_size * n
    logging.info(f"Splitting genome into {n} files, each with {split_size} bases")
    logging.info(f"Total genome length: {total_len} bases")
    logging.info(f"Total bases covered: {bases_will_be_covered} bases")
    logging.info(f"Final split will be {bases_will_be_covered - total_len} bases shorter than the other splits.")
    
    cur_index = 0
    cur_split_covered = 0
    for chrom, length in chroms:
        cur_start = 0
        while True:
            bp_left_in_chrom = length - cur_start
            bp_left_in_split = split_size - cur_split_covered
            if bp_left_in_chrom < bp_left_in_split:
                cur_end = length
            else:
                cur_end = cur_start + bp_left_in_split
            # write to file
            outfiles[cur_index].write(f"{chrom}\t{cur_start}\t{cur_end}\n")

            # reset for next iteration
            cur_split_covered += cur_end - cur_start
            cur_start = cur_end

            # we have filled up the split
            if cur_split_covered == split_size:
                outfiles[cur_index].close()
                cur_split_covered = 0
                cur_index += 1

            # we have used up the entire chromosome
            if cur_start == length:
                break

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

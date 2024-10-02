#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional
import pysam
import pandas as pd
import tqdm


def main(
    vcf: Optional[Path] = None,
    tbl: Optional[Path] = None,
    o_vcf: Optional[Path] = None,
    *,
    verbose: int = 3,
    threads: int = 8,
    max_p_values: float = 5e-8,
    max_1kg_af: float = 1.0,
):
    """
    Author Mitchell R. Vollger

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    :param count: Number of times to display the greeting
    :param verbose: Set the logging level of the function
    """
    if o_vcf is None:
        o_vcf = sys.stdout

    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    gwas = pd.read_csv(tbl, sep="\t")
    logging.info(f"Read in {len(gwas):,} GWAS SNPs")
    # filter for rows where the rsID starts with rs
    # filter out rsIDs that don't end with a number
    gwas = gwas[
        gwas["rsID"].str.startswith("rs") & gwas["rsID"].str[2:].str.isnumeric()
    ]
    logging.info(f"Filtered to {len(gwas):,} for formatting of rsIDs")
    gwas = gwas[(gwas["P_VAL"] < max_p_values)]
    logging.info(f"Filtered to {len(gwas):,} GWAS SNPs with p-value < {max_p_values}")
    # strip the leading rs from the rsID and convert to an integer
    gwas["rsID"] = gwas["rsID"].str.lstrip("rs").astype(int)
    # filter for unique rsIDs
    gwas = gwas.drop_duplicates("rsID")
    # make the index the rsID
    gwas = gwas.set_index("rsID")
    logging.info(f"Filtered to {len(gwas):,} unique GWAS SNPs (rsID and risk allele)")

    # read in the VCF file with pysam
    vcf = pysam.VariantFile(vcf, threads=threads)
    o_vcf = pysam.VariantFile(o_vcf, "wb", header=vcf.header, threads=threads)
    logging.info(f"Reading in VCF file {vcf}")
    count = 0
    # loop over the VCF file
    for rec in tqdm.tqdm(vcf.fetch()):
        # get the rsID from the INFO field
        rsID = rec.info.get("rsID", -1)
        if rsID != -1 and rsID in gwas.index:
            risk_alt = gwas.loc[rsID, "ALT"]
            alt = rec.alts[0]
            if risk_alt == alt:
                # get the 1KG allele frequency
                af = float(rec.info.get("FREQ", ("1000Genomes:1.0"))[0].split(":")[1])
                if af > max_1kg_af:
                    continue
                count += 1
                o_vcf.write(rec)
    o_vcf.close()
    logging.info(f"Found {count:,} GWAS risk alleles in VCF")
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")

#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional


def main(
    *,
    endpoint: str = "kopah",
    verbose: int = 3,
):
    """
    Convert S3 locations from stdin to URLs on stdout
    
    Author Mitchell R. Vollger

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    :param endpoint: S3 endpoint format (kopah, aws, or custom URL)
    :param verbose: Set the logging level of the function
    """
    infile=sys.stdin
    outfile=sys.stdout  
    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    # Define endpoint formats
    def get_url(s3_path):
        path = s3_path[5:]  # Remove s3:// prefix
        
        if endpoint == "kopah":
            return f"https://s3.kopah.uw.edu/{path}"
        elif endpoint == "aws":
            # Split into bucket and key for AWS format
            parts = path.split('/', 1)
            if len(parts) >= 2:
                bucket, key = parts[0], parts[1]
                return f"https://{bucket}.s3.amazonaws.com/{key}"
            else:
                bucket = parts[0]
                return f"https://{bucket}.s3.amazonaws.com/"
        else:
            # Custom endpoint - assume it's a base URL
            if not endpoint.startswith('http'):
                endpoint_url = f"https://{endpoint}"
            else:
                endpoint_url = endpoint
            return f"{endpoint_url}/{path}"

    count = 0
    for line in infile:
        s3_location = line.strip()
        if not s3_location:
            continue
            
        if s3_location.startswith('s3://'):
            url = get_url(s3_location)
            outfile.write(url + '\n')
            count += 1
        else:
            logging.warning(f"Skipping non-S3 location: {s3_location}")
    

    logging.info(f"Converted {count:,} S3 locations to URLs")
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
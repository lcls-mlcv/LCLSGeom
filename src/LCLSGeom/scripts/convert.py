"""
Command line tool to convert files between different formats.
"""

import sys
import argparse
import logging

logger = logging.getLogger(__name__)

from LCLSGeom.converter import PsanaToCrystFEL, CrystFELToPsana

def main():
    parser = argparse.ArgumentParser(description="Convert geometry files between different formats.")

    # Required arguments
    parser.add_argument("-i", "--in_file", type=str, help="Path to the input file")
    parser.add_argument("-d", "--detname", type=str, help="Detector name")
    parser.add_argument("-o", "--out_file", type=str, help="Path to the output file")
    args = parser.parse_args()

    if args.in_file.endswith(".data") and args.out_file.endswith(".geom"):
        PsanaToCrystFEL.convert(in_file=args.in_file, detname=args.detname, out_file=args.out_file)
    elif args.in_file.endswith(".geom") and args.out_file.endswith(".data"):
        CrystFELToPsana.convert(in_file=args.in_file, detname=args.detname, out_file=args.out_file)
    else:
        logger.error("Unsupported conversion between %s and %s.", args.in_file, args.out_file)
        sys.exit(1)
    logger.info("Geometry converted successfully.")

if __name__ == "__main__":
    main()
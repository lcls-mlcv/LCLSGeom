"""
Command line tool to push calibration constants to the calibration database. 
This is a wrapper around the psana tools, and is meant to be used in a containerized environment where psana is available.
"""

import sys
import argparse
import logging

logger = logging.getLogger(__name__)

import psana

if hasattr(psana, "xtc_version"):
    IS_PSANA2 = True
else:
    IS_PSANA2 = False
    logger.error("This is a psana2 tool. Please source the appropriate environment.")
    sys.exit(1)

from LCLSGeom.manager import push_to_database

def main():
    parser = argparse.ArgumentParser(description="Push geometry file to calibration database.")
    parser.add_argument("-e", "--exp", type=str, help="Experiment name")
    parser.add_argument("-r", "--run", type=int, help="Run number")
    parser.add_argument("-d", "--detname", type=str, help="Detector name")
    parser.add_argument("-g", "--geometry", type=str, help="Path to the geometry file to push")
    args = parser.parse_args()

    push_to_database(exp=args.exp, run=args.run, detname=args.detname, out_file=args.geometry)
    logger.info("Geometry pushed to database successfully.")

if __name__ == "__main__":
    main()

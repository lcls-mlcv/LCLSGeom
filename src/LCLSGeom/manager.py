"""
LCLS-II detectors:
File manager functions to handle push and pull of geometry files to and from the psana calibration database.
This module also includes functions to fetch default geometry files if not found in the calibration database.

LCLS-I detectors:
File manager functions to look for geometry in the experiment calibration directory.
This module also includes functions to fetch default geometry files if not found in the experiment calibration directory.
"""

import logging

import psana
if hasattr(psana, "xtc_version"):
    from psana import DataSource
    import psana.pscalib.calib.MDBUtils as mu
    import psana.pscalib.calib.MDBWebUtils as wu
    import psana.detector.UtilsCalib as uc
    cc = wu.cc

    IS_PSANA2 = True
else:
    IS_PSANA2 = False

import os
import re
import tempfile
from typing import Optional

from LCLSGeom.utils.calibmap import group_from_detname, source_from_detname, select_calib_file, get_full_detname

from logging import getLogger

logger: logging.Logger = getLogger(__name__)


def get_geometry(detname: str, exp: Optional[str]=None, run: Optional[int]=None) -> str:
    """
    Get the geometry file for a given experiment, run, and detector name. 
    If LCLS-II, check first the calibration database, then the default geometry file.
    If LCLS-I, check first the experiment calibration directory, then the default geometry file.
    If not found in either location, get default geometry file from templates.

    Parameters
    ----------
    exp : str
        Experiment name
    run : int
        Run number
    detname : str
        Name of the detector

    Returns
    -------
    in_file : str
        Path to the geometry file
    """
    if exp and run:
        if IS_PSANA2:
            in_file = pull_from_database(exp=exp, run=run, detname=detname)
            if in_file:
                logger.info(f"Geometry file for {detname} found in calibration database for experiment {exp} and run {run}.")
        else:
            in_file = check_calibration_directory(exp=exp, run=run, detname=detname)
            if detname.lower() == "rayonix":
                in_file = get_default_Rayonix_geometry(exp=exp, run=run)
            if in_file:
                logger.info(f"Geometry file for {detname} found in experiment calibration directory for experiment {exp} and run {run}.")
    
        if not in_file:
            logger.info(f"Fetching default geometry file for {detname}.")
            in_file = get_default_geometry(detname=detname)
    else:
        logger.info(f"Fetching default geometry file for {detname}.")
        in_file = get_default_geometry(detname=detname)
    return in_file


def get_default_geometry(detname: str) -> str:
    """
    Get the default geometry file for a given detector name.
    For LCLS-II, check first the psana installation for the default geometry file, then fetch from templates if not found.

    Parameters
    ----------
    detname : str
        Name of the detector in psana

    Returns
    -------
    in_file : str
        Path to the default geometry file
    """
    detname = get_full_detname(detname)
    if IS_PSANA2:
        psana_path = psana.__path__[0]
        data_dir = os.path.join(psana_path, "pscalib", "geometry", "data")
        in_file = os.path.join(data_dir, f"geometry-def-{detname}.data")
        if not os.path.exists(in_file):
            logger.warning(f"Default geometry file not found for {detname} in psana installation. Fetching from templates.")
            in_file = fetch_template(detname=detname)

    else:
        in_file = fetch_template(detname=detname)

    if not os.path.exists(in_file):
        raise FileNotFoundError(f"Default geometry file not found for detector {detname}.")

    return in_file


def pull_from_database(exp: str, run: int, detname: str) -> Optional[str]:
    """
    For LCLS-II, pull the geometry file for a given experiment, run, and detector name from the calibration database.
    
    Parameters
    ----------
    exp : str
        Experiment name
    run : int
        Run number
    detname : str
        Name of the detector

    Returns
    -------
    in_file : str
        Path to the geometry file if found written to scratch space, else None
    """
    ds = DataSource(exp=exp, run=run)
    runs = next(ds.runs())
    detector = runs.Detector(detname)
    longname: str = detector.raw._uniqueid
    shortname: str = uc.detector_name_short(longname)
    results = wu.calib_constants(
        shortname,
        exp=exp,
        run=run,
        ctype="geometry",
        url=cc.URL,
    )
    if results is None:
        return None

    data, _ = results
    with tempfile.NamedTemporaryFile(mode="w", suffix=".data", delete=False) as in_file:
        in_file.write(data)
        return in_file.name


def push_to_database(exp: str, run: int, detname: str, out_file: str) -> None:
    """
    Upload the geometry to the calibration database.

    Parameters
    ----------
    exp : str
        Experiment name
    run : int
        Run number
    detname : str
        Name of the detector
    out_file : str
        Path to the geometry file to be uploaded
    """
    ctype = "geometry"
    dtype = "str"
    data = mu.data_from_file(out_file, ctype, dtype)
    ds = DataSource(exp=exp, run=run, skip_calib_load="all", max_events=1)
    runs = next(ds.runs())
    detector = runs.Detector(detname)
    longname: str = detector.raw._uniqueid
    shortname: str = uc.detector_name_short(longname)
    det_type: str = detector._dettype
    run_orig: int = run
    run_beg: int = run
    run_end: str = "end"
    kwa = {
        "iofname": out_file,
        "experiment": exp,
        "ctype": ctype,
        "dtype": dtype,
        "detector": shortname,
        "shortname": shortname,
        "detname": detname,
        "longname": longname,
        "run": run,
        "run_beg": run_beg,
        "run_end": run_end,
        "run_orig": run_orig,
        "dettype": det_type,
    }
    _ = wu.deploy_constants(
        data,
        exp,
        longname,
        url=cc.URL_KRB,
        krbheaders=cc.KRBHEADERS,
        **kwa,
    )


def check_calibration_directory(exp: str, run: int, detname: str) -> Optional[str]:
    """
    Check for the geometry file for a given experiment, run, and detector name in the experiment calibration directory.

    Parameters
    ----------
    exp : str
        Experiment name
    run : int
        Run number
    detname : str
        Name of the detector
    
    Returns
    -------
    in_file : str
        Path to the geometry file if found, else None
    """
    cdir = f"/sdf/data/lcls/ds/{exp[:3]}/{exp}/calib"
    src = source_from_detname(detname, exp[:3])
    group = group_from_detname(detname)
    ctype = "geometry"
    calib_dir = os.path.join(cdir, src, group, ctype)
    in_file = select_calib_file(calib_dir, run)
    return in_file


def fetch_template(detname):
    """
    Get the template geometry file for a given detector name

    Parameters
    ----------
    detname : str
        Detector name

    Returns
    -------
    in_file : str
        Path to the template geometry file
    """
    current_dir = os.path.dirname(__file__)
    in_file = os.path.join(current_dir, "templates", f"geometry-def-{detname}.data")

    if not os.path.exists(in_file):
        raise FileNotFoundError(f"Detector {detname} not supported.")

    return in_file


def update_Rayonix_binning(in_file: str, binning_fast: int, binning_slow: int) -> str:
    """
    Update binning in a Rayonix MTRX:V2 geometry file with the corresponding shape and pixel size.

    Parameters
    ----------
    file_path : str
        Path to the geometry file to update
    binning_fast : int
        Fast binning factor
    binning_slow : int
        Slow binning factor
        Pixel size in micrometers
    """
    shape_slow = int(7680 / binning_slow)
    shape_fast = int(7680 / binning_fast)
    pixel_size_slow = int(44 * binning_slow)
    pixel_size_fast = int(44 * binning_fast)

    with open(in_file, "r") as f:
        lines = f.readlines()

    pattern = re.compile(r"(MTRX:V2:)\d+:\d+:\d+:\d+")
    replacement = f"MTRX:V2:{shape_slow}:{shape_fast}:{pixel_size_slow}:{pixel_size_fast}"

    for i, line in enumerate(lines):
        if "MTRX:V2" in line:
            lines[i] = pattern.sub(replacement, line)
            break
    
    data = "".join(lines)
    return data 


def get_default_Rayonix_geometry(exp: str, run: int) -> str:
    """
    Get the default geometry file for the LCLS-I Rayonix detector with the appropriate binning based on pixel size.

    Parameters
    ----------
    exp : str
        Experiment name
    run : int
        Run number


    Returns
    -------
    in_file : str
        Path to the default geometry file with updated binning for the Rayonix detector
    """
    ds = psana.DataSource(f"exp={exp}:run={run}:idx")
    env = ds.env()
    cfg = env.configStore()
    binning_fast = cfg.get(psana.Rayonix.ConfigV2).binning_f()
    binning_slow = cfg.get(psana.Rayonix.ConfigV2).binning_s()
    in_file = get_default_geometry(detname="rayonix")
    data = update_Rayonix_binning(in_file=in_file, binning_fast=binning_fast, binning_slow=binning_slow)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".data", delete=False) as tmp_file:
        tmp_file.write(data)
        in_file = tmp_file.name
    return in_file
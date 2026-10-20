"""
File manager functions to locate a geometry file for a detector.

For LCLS-II the geometry is pushed to and pulled from the psana calibration database; for
LCLS-I it is looked up in the experiment calibration directory. Either way, a shipped default
template is used if nothing is found.
"""

import psana

IS_PSANA2 = hasattr(psana, "xtc_version")
if IS_PSANA2:
    from psana import DataSource
    import psana.pscalib.calib.MDBUtils as mu
    import psana.pscalib.calib.MDBWebUtils as wu
    import psana.detector.UtilsCalib as uc
    cc = wu.cc

import logging
import os
import re
import tempfile
from typing import Optional

from LCLSGeom.utils.pscalib_utils import (get_full_detname, group_from_detname,
                                          select_calib_file, source_from_detname)

logger: logging.Logger = logging.getLogger(__name__)


def get_geometry(detname: str, exp: Optional[str]=None, run: Optional[int]=None) -> str:
    """
    Get the geometry file for a given detector, experiment and run.

    If LCLS-II, check first the calibration database, then the default geometry file.
    If LCLS-I, check first the experiment calibration directory, then the default geometry file.
    Without an experiment and run, the default geometry file is used directly.

    The LCLS-I Rayonix is a special case: its geometry depends on the binning configured for the
    run, so it is always built from the template rather than read from the calibration directory.

    Parameters
    ----------
    detname : str
        Name of the detector
    exp : str
        Experiment name
    run : int
        Run number

    Returns
    -------
    in_file : str
        Path to the geometry file
    """
    in_file = None
    if exp and run:
        if IS_PSANA2:
            in_file = pull_from_database(exp=exp, run=run, detname=detname)
            source = "calibration database"
        elif detname.lower() == "rayonix":
            in_file = get_default_Rayonix_geometry(exp=exp, run=run)
            source = "Rayonix template with run binning"
        else:
            in_file = check_calibration_directory(exp=exp, run=run, detname=detname)
            source = "experiment calibration directory"
        if in_file:
            logger.info(f"Geometry file for {detname} found in {source} for experiment {exp} and run {run}.")

    if not in_file:
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
        Path to a temporary file holding the geometry, or None if the database has none
    """
    ds = DataSource(exp=exp, run=run)
    detector = next(ds.runs()).Detector(detname)
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
    detector = next(ds.runs()).Detector(detname)
    longname: str = detector.raw._uniqueid
    shortname: str = uc.detector_name_short(longname)
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
        "run_beg": run,
        # The geometry stays valid from this run on, until a later deployment supersedes it
        "run_end": "end",
        "run_orig": run,
        "dettype": detector._dettype,
    }
    wu.deploy_constants(
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
    hutch = exp[:3]
    calib_dir = os.path.join(
        f"/sdf/data/lcls/ds/{hutch}/{exp}/calib",
        group_from_detname(detname),
        source_from_detname(detname, hutch),
        "geometry",
    )
    return select_calib_file(calib_dir, run)


def fetch_template(detname):
    """
    Get the shipped template geometry file for a given detector name

    Parameters
    ----------
    detname : str
        Full detector name, as returned by get_full_detname

    Returns
    -------
    in_file : str
        Path to the template geometry file
    """
    in_file = os.path.join(os.path.dirname(__file__), "templates", f"geometry-def-{detname}.data")

    if not os.path.exists(in_file):
        raise FileNotFoundError(f"Detector {detname} not supported.")

    return in_file


def update_Rayonix_binning(in_file: str, binning_fast: int, binning_slow: int) -> str:
    """
    Update binning in a Rayonix MTRX:V2 geometry file with the corresponding shape and pixel size.

    Parameters
    ----------
    in_file : str
        Path to the geometry file to update
    binning_fast : int
        Fast binning factor
    binning_slow : int
        Slow binning factor

    Returns
    -------
    data : str
        Contents of the geometry file with the binning updated
    """
    shape_slow = int(7680 / binning_slow)
    shape_fast = int(7680 / binning_fast)
    pixel_size_slow = int(44 * binning_slow)
    pixel_size_fast = int(44 * binning_fast)

    pattern = re.compile(r"MTRX:V2:\d+:\d+:\d+:\d+")
    replacement = f"MTRX:V2:{shape_slow}:{shape_fast}:{pixel_size_slow}:{pixel_size_fast}"
    with open(in_file, "r") as f:
        return pattern.sub(replacement, f.read(), count=1)


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
        Path to a temporary geometry file with the binning of the run
    """
    ds = psana.DataSource(f"exp={exp}:run={run}:idx")
    cfg = ds.env().configStore().get(psana.Rayonix.ConfigV2)
    data = update_Rayonix_binning(
        in_file=get_default_geometry(detname="rayonix"),
        binning_fast=cfg.binning_f(),
        binning_slow=cfg.binning_s(),
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".data", delete=False) as tmp_file:
        tmp_file.write(data)
        return tmp_file.name
"""
PSCalib utility functions and mappings for psana detector names.
"""

import os
from math import atan2, degrees
from time import strftime
from typing import Optional
import getpass

import numpy as np

# psana detname -> (template stem, calibration source).
# One row per detector rather than index-coupled parallel tuples: a row is read as a whole, so a
# missing entry can no longer shift one column out of alignment with the others. The calibration
# group is not stored because it follows from the source; see group_from_detname.
PSANA_DETECTORS = {
    "rayonix":       ("rayonix",       "Rayonix"),
    "jungfrau05m":   ("jungfrau05M",   "Jungfrau"),
    "jungfrau1m":    ("jungfrau1M",    "Jungfrau"),
    "jungfrau4m":    ("jungfrau4M",    "Jungfrau"),
    # psana calls the 16M plain "jungfrau"
    "jungfrau":      ("jungfrau16M",   "Jungfrau"),
    "epix10k2m":     ("epix10k2M",     "Epix10ka2M"),
    "epix10kaquad0": ("epix10kaQuad0", "Epix10kaQuad.0"),
    "epix10kaquad1": ("epix10kaQuad1", "Epix10kaQuad.1"),
    "epix10kaquad2": ("epix10kaQuad2", "Epix10kaQuad.2"),
    "epix10kaquad3": ("epix10kaQuad3", "Epix10kaQuad.3"),
}

# Calibration group per detector type, keyed by the source name with any quad index stripped
CALIB_GROUPS = {
    "Rayonix": "Camera::CalibV1",
    "Jungfrau": "Jungfrau::CalibV1",
    "Epix10ka2M": "Epix10ka2M::CalibV1",
    "Epix10kaQuad": "Epix10kaQuad::CalibV1",
}

hutch_to_station = {
    "XPP": "XppEndstation.0",
    "XCS": "XcsEndstation.0",
    "CXI": "CxiDs1.0",
    "MEC": "MecTargetChamber.0",
    "MFX": "MfxEndstation.0",
}


def _lookup(detname: str) -> tuple:
    """Row of PSANA_DETECTORS for a psana detector name."""
    row = PSANA_DETECTORS.get(detname.lower())
    if row is None:
        raise ValueError(f"Unknown detector type: {detname}")
    return row


def get_full_detname(detname: str) -> str:
    """Start from psana detector name and return the full detector name."""
    return _lookup(detname)[0]


def get_calib_detname(detname: str) -> str:
    """Start from psana detector name and return the calibration detector name."""
    return _lookup(detname)[1]


def group_from_detname(detname: str) -> str:
    """Retrieve the calibration group string from the detector name."""
    # The quad index is part of the source name but not of the group name
    return CALIB_GROUPS.get(get_calib_detname(detname).split(".")[0])


def source_from_detname(detname: str, hutch: str) -> str:
    """Retrieve the source string from the detector name and hutch."""
    station = hutch_to_station.get(hutch.upper())
    if station is None:
        raise ValueError(f"Unknown hutch: {hutch}")
    calib_detname = get_calib_detname(detname)
    # The quad index is already part of the calibration name, so no trailing .0 is appended
    if "Epix10kaQuad" in calib_detname:
        return f"{station}:{calib_detname}"
    return f"{station}:{calib_detname}.0"


def select_calib_file(calib_dir: str, run: int) -> Optional[str]:
    """
    Select the calibration file covering a run in a calibration directory.

    Calibration files are named '<run-begin>-<run-end>.data', where <run-end> is either a run
    number or 'end'. The highest-numbered range covering the run wins; if none does, the
    catch-all '0-end.data' is returned whether or not it exists.

    Parameters
    ----------
    calib_dir : str
        Path to the calibration directory
    run : int
        Run number

    Returns
    -------
    in_file : str
        Path to the calibration file, or None if the directory does not exist
    """
    if not os.path.exists(calib_dir):
        return

    RUN_MAX = 9999
    ranges = []
    for fname in os.listdir(calib_dir):
        basename, ext = os.path.splitext(fname)
        if ext != ".data":
            continue
        fields = basename.split("-")
        if len(fields) != 2:
            continue
        begin, end = fields
        if not begin.isdigit() or not (end.isdigit() or end == "end"):
            continue
        end_int = RUN_MAX if end == "end" else int(end)
        ranges.append((int(begin), end_int, os.path.join(calib_dir, fname)))

    for begin, end, file in sorted(ranges, reverse=True):
        if begin <= run <= end:
            return file

    return os.path.join(calib_dir, "0-end.data")


def shape3d(geo):
    """
    Return 3d shape of the arrays as (<number-of-segments>, <rows>, <cols>)

    Parameters
    ----------
    geo : GeometryAccess
        psana geometry object

    Returns
    -------
    shape : tuple
        Calibration shape as (<number-of-segments>, <rows>, <cols>)
    """
    seg = geo.get_seg_geo().algo
    sshape = seg.shape()
    x, _, _ = geo.get_pixel_coords()
    return (int(x.size/seg.size()), sshape[0], sshape[1])

def get_root_geometry(geo, image_frame=True):
    """
    Return the camera node of a pscalib geometry, optionally put into the image frame

    In the image frame the camera node offsets are zeroed and the 90 degree rotation of the
    psana convention is undone, leaving:
    - X-axis horizontal from left to right,
    - Y-axis vertical from bottom to top,
    - Z-axis opposite the beam direction.

    Parameters
    ----------
    geo : GeometryAccess
        psana geometry object, as returned by pscalib GeometryAccess
    image_frame : bool
        If True, zero the camera node placement; otherwise, leave the psana laboratory frame

    Returns
    -------
    camera : GeometryObject
        The camera node, with its placement zeroed if image_frame is True
    """
    camera = geo.get_top_geo().get_list_of_children()[0]
    if image_frame:
        camera.x0 = 0
        camera.y0 = 0
        camera.z0 = 0
        camera.rot_x = 0
        camera.rot_y = 0
        camera.rot_z -= 90
    return camera


def angle_and_tilt(a):
    """
    Split an angle into the nearest psana design angle and the residual tilt

    psana stores a segment orientation as a design rotation, which is a multiple of 90 degrees,
    plus a small tilt. The returned angle is wrapped into [0, 360).

    Parameters
    ----------
    a : float
        Angle in degrees

    Returns
    -------
    angle, tilt : float, float
        Design angle and residual tilt, in degrees
    """
    desangles = np.array((-180,-90, 0, 90, 180))
    difangles = a-desangles
    absdifang = np.abs(difangles)
    imin = np.argmin(absdifang)
    angle, tilt = desangles[imin], difangles[imin]
    return (angle if angle>=0 else angle+360), tilt

def unit_vector_pitch_angle_max_ind(u):
    """
    Pitch angle of a unit vector out of the xy plane, and the axis it lies closest to

    The angle is wrapped into [-90, 90], so the result is independent of the vector's sign.

    Parameters
    ----------
    u : numpy.ndarray
        Input unit vector

    Returns
    -------
    pitch, imax : float, int
        Pitch angle in degrees, and the index of the largest component of the vector
    """
    absu = np.abs(u)
    imax = np.where(absu == np.amax(absu))[0]
    pitch = degrees(atan2(u[2],u[imax]))
    pitch = (pitch+180) if pitch<-90 else (pitch-180) if pitch>90 else pitch
    return pitch, imax

def tilt_xy(uf, us):
    """
    Tilts about the x and y axes of the plane spanned by the fast and slow scan vectors

    Which vector carries which tilt depends on the axis the fast scan vector runs along.

    Parameters
    ----------
    uf : numpy.ndarray
        Unit vector along the fast scan direction
    us : numpy.ndarray
        Unit vector along the slow scan direction

    Returns
    -------
    tilt_x, tilt_y : float, float
        Tilts about the x and y axes, in degrees
    """
    tilt_f, imaxf = unit_vector_pitch_angle_max_ind(uf)
    tilt_s, _ = unit_vector_pitch_angle_max_ind(us)
    tilt_x, tilt_y = (tilt_s, tilt_f) if imaxf == 0 else (tilt_f, tilt_s)
    return tilt_x, -tilt_y


# title, DETECTOR field, and COMMENT subject, per detector. '{detname}' is the name as passed.
PSANA_HEADERS = {
    'rayonix':      ('Rayonix',      'Rayonix',                            'the Rayonix'),
    'epix10k2m':    ('ePix10k2M',    'Epix10ka2M',                         'the 16-segment ePix10k2M'),
    'epix10kaquad': ('ePix10kaQuad', '{detname}',                          'the ePix10kaQuad'),
    'jungfrau05m':  ('Jungfrau1M',   'DetLab.0:Jungfrau.0 or jungfrau05M', 'the 1-segment Jungfrau1M'),
    'jungfrau1m':   ('Jungfrau1M',   'DetLab.0:Jungfrau.1 or jungfrau1M',  'the 2-segment Jungfrau1M'),
    'jungfrau4m':   ('Jungfrau4M',   'DetLab.0:Jungfrau.2 or jungfrau4M',  'the 8-segment Jungfrau4M'),
    'jungfrau16m':  ('Jungfrau16M',  'DetLab.0:Jungfrau.3 or jungfrau16M', 'the 32-segment Jungfrau16M'),
}

def header_psana(detname):
    """
    Returns the header for the psana geometry file based on the detector name.
    The header includes information such as the title, date, author, experiment,
    detector, calibration type, and comments.

    An unrecognized detector name yields the parameter block alone, with no title.

    Parameters
    ----------
    detname : str
        The name of the detector (e.g., 'Rayonix', 'ePix10k2M', etc.)
    """
    entry = PSANA_HEADERS.get(detname.lower())
    if entry is None:
        txt = ''
    else:
        title, detector, subject = entry
        txt = (
            '# TITLE      Geometry parameters of %s' % title
            + '\n# DATE_TIME  %s' % strftime('%Y-%m-%dT%H:%M:%S %Z')
            + '\n# METROLOGY  no metrology available'
            + '\n# AUTHOR     %s' % getpass.getuser()
            + '\n# EXPERIMENT N/A'
            + '\n# DETECTOR   %s' % detector.format(detname=detname)
            + '\n# CALIB_TYPE geometry'
            + '\n# COMMENT:01 Automatically created from BayFAI for %s detector' % subject
            + '\n'
        )
    txt +=\
        '\n# PARAM:01 PARENT     - name and version of the parent object'\
        '\n# PARAM:02 PARENT_IND - index of the parent object'\
        '\n# PARAM:03 OBJECT     - name and version of the object'\
        '\n# PARAM:04 OBJECT_IND - index of the new object'\
        '\n# PARAM:05 X0         - x-coordinate [um] of the object origin in the parent frame'\
        '\n# PARAM:06 Y0         - y-coordinate [um] of the object origin in the parent frame'\
        '\n# PARAM:07 Z0         - z-coordinate [um] of the object origin in the parent frame'\
        '\n# PARAM:08 ROT_Z      - object design rotation angle [deg] around Z axis of the parent frame'\
        '\n# PARAM:09 ROT_Y      - object design rotation angle [deg] around Y axis of the parent frame'\
        '\n# PARAM:10 ROT_X      - object design rotation angle [deg] around X axis of the parent frame'\
        '\n# PARAM:11 TILT_Z     - object tilt angle [deg] around Z axis of the parent frame'\
        '\n# PARAM:12 TILT_Y     - object tilt angle [deg] around Y axis of the parent frame'\
        '\n# PARAM:13 TILT_X     - object tilt angle [deg] around X axis of the parent frame'\
        '\n'\
        '\n# HDR PARENT IND     OBJECT IND    X0[um]   Y0[um]   Z0[um]   ROT-Z  ROT-Y  ROT-X     TILT-Z    TILT-Y    TILT-X'
    return txt

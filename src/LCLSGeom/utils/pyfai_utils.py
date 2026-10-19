"""
pyFAI utility functions for loading a .poni file, applying its PONI parameters, and moving
coordinates between the image/psana-lab frames and the pyFAI frame.
"""

import numpy as np

from LCLSGeom.utils.detectors import get_family


def get_pyFAI_detname(detname):
    """
    Start from a psana detector name and return the pyFAI detector name.

    pyFAI names a detector by its family, so this is the family, unchanged.

    Parameters
    ----------
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)

    Returns
    -------
    pyFAI_detname : str
        Detector name in pyFAI (e.g., 'Jungfrau4M', 'ePix10k2M', etc.)
    """
    return get_family(detname)


def load_poni(poni_file):
    """
    Load a PyFAI .poni file and return its geometry parameters

    PyFAI is imported lazily so that the psana, CrystFEL and dxtbx converters do not require it.

    Parameters
    ----------
    poni_file : str
        Path to the PyFAI .poni file

    Returns
    -------
    params : list
        PONI parameters as (distance, poni1, poni2, rot1, rot2, rot3)
    """
    import pyFAI

    pyFAI.use_opencl = False
    return pyFAI.load(poni_file).param


def pyfai_to_image_poni(detector, in_file, out_file):
    """
    Re-express a .poni fitted on a PsanaToPyFAI detector in the det.image frame (pyFAI orientation 2)

    Parameters
    ----------
    detector : pyFAI.Detector
        PsanaToPyFAI detector the .poni was fitted on, built with image_frame=True
    in_file : str
        Path to the input .poni file, relative to the metrology
    out_file : str
        Path to the output .poni file, relative to the assembled image
    """
    import pyFAI
    from pyFAI.detectors import Detector
    from pyFAI.geometry import Geometry

    ai = pyFAI.load(in_file)
    rows, cols = detector.pixel_index_map[..., 0], detector.pixel_index_map[..., 1]
    p1, p2, _ = detector.calc_cartesian_positions()
    pix = detector.pixel1
    # The median absorbs the integer snapping of pixel_index_map
    y0 = np.median(p1.reshape(rows.shape) + (rows + 0.5) * pix)
    x0 = np.median(p2.reshape(cols.shape) - (cols + 0.5) * pix)
    shape = (int(rows.max()) + 1, int(cols.max()) + 1)
    Geometry(
        dist=ai.dist, poni1=ai.poni1 - (y0 - shape[0] * pix), poni2=ai.poni2 - x0,
        rot1=ai.rot1, rot2=ai.rot2, rot3=ai.rot3, wavelength=ai.wavelength,
        detector=Detector(pix, pix, max_shape=shape, orientation=2),
    ).save(out_file)


def rotation_matrix(params):
    """
    Compose the three PONI rotations into a single rotation matrix

    Parameters
    ----------
    params : list
        Detector parameters found by PyFAI calibration, rot1/rot2/rot3 at indices 3, 4 and 5

    Returns
    -------
    matrix : np.ndarray
        3x3 rotation matrix taking detector coordinates to sample coordinates
    """
    cos_rot1 = np.cos(params[3])
    cos_rot2 = np.cos(params[4])
    cos_rot3 = np.cos(params[5])
    sin_rot1 = np.sin(params[3])
    sin_rot2 = np.sin(params[4])
    sin_rot3 = np.sin(params[5])
    # Rotation about axis 1: note this rotation is left-handed
    rot1 = np.array([[1.0, 0.0, 0.0],
                     [0.0, cos_rot1, sin_rot1],
                     [0.0, -sin_rot1, cos_rot1]])
    # Rotation about axis 2: note this rotation is left-handed
    rot2 = np.array([[cos_rot2, 0.0, -sin_rot2],
                     [0.0, 1.0, 0.0],
                     [sin_rot2, 0.0, cos_rot2]])
    # Rotation about axis 3: note this rotation is right-handed
    rot3 = np.array([[cos_rot3, -sin_rot3, 0.0],
                     [sin_rot3, cos_rot3, 0.0],
                     [0.0, 0.0, 1.0]])
    return rot3 @ rot2 @ rot1


def correct_geometry(detector, params):
    """
    Apply the PONI parameters to the detector pixel coordinates

    The PONI translation is subtracted from the pixel coordinates and the PONI rotations are
    applied, bringing the pixels from the detector frame into the sample frame.

    Parameters
    ----------
    detector : pyFAI.detectors.Detector
        The pyFAI Detector object to be corrected
    params : list
        Detector parameters found by PyFAI calibration

    Returns
    -------
    x, y, z : np.ndarray
        Flat arrays of corrected pixel coordinates in meters, in the sample frame
    """
    p1, p2, p3 = detector.calc_cartesian_positions()
    dist = params[0]
    poni1 = params[1]
    poni2 = params[2]
    p1 = (p1 - poni1).ravel()
    p2 = (p2 - poni2).ravel()
    if p3 is None:
        p3 = np.zeros_like(p1) + dist
    else:
        p3 = (p3+dist).ravel()
    coord_det = np.stack((p1, p2, p3), axis=0)
    coord_sample = np.dot(rotation_matrix(params), coord_det)
    x, y, z = coord_sample
    return x, y, z


def image_to_pyfai(x, y, z, image_frame=True):
    """
    Convert image or psana-lab coordinates to pyFAI coordinates

    z is recentered on its mean, since pyFAI carries the sample-detector distance separately,
    and is flattened to zero for a detector whose segments all sit at the same z.

    Parameters
    ----------
    x : np.ndarray
        X coordinate in microns
    y : np.ndarray
        Y coordinate in microns
    z : np.ndarray
        Z coordinate in microns
    image_frame : bool
        If True, x, y, z are image frame coordinates; otherwise, x, y, z are psana laboratory frame coordinates

    Returns
    -------
    x, y, z : np.ndarray
        pyFAI coordinates in meters
    """
    x = x * 1e-6
    y = y * 1e-6
    z = z * 1e-6
    if len(np.unique(z))==1:
        z = np.zeros_like(z)
    else:
        z -= np.mean(z)

    if image_frame:
        return y, x, -z
    else:
        return -x, y, -z

def pyfai_to_image(x, y, z, image_frame=True):
    """
    Convert pyFAI coordinates back to image or psana-lab coordinates

    Parameters
    ----------
    x : np.ndarray
        X coordinate in meters
    y : np.ndarray
        Y coordinate in meters
    z : np.ndarray
        Z coordinate in meters
    image_frame : bool
        If True, return image frame coordinates; otherwise, return psana laboratory frame coordinates

    Returns
    -------
    x, y, z : np.ndarray
        Image or psana laboratory frame coordinates in microns
    """
    x = x * 1e6
    y = y * 1e6
    z = z * 1e6

    if image_frame:
        return y, x, -z
    else:
        return -x, y, -z

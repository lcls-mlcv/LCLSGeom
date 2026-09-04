"""
Frame transformations between image/psana-lab coordinates and pyFAI coordinates.
"""

import numpy as np

def image_to_pyfai(x, y, z, image_frame=True):
    """
    Convert image coordinates to pyfai coordinates

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
    Convert back to image coordinates

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

    """
    x = x * 1e6
    y = y * 1e6
    z = z * 1e6

    if image_frame:
        return y, x, -z
    else:
        return -x, y, -z
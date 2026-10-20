"""
Geometry helpers shared by every converter.
"""

import numpy as np

def asic_axes(x, y, z, r0, c0, arows, acols):
    """
    Unit vectors along the fast and slow scan directions of one ASIC

    Parameters
    ----------
    x, y, z : np.ndarray
        Pixel coordinates of a single segment, shaped (<rows>, <cols>)
    r0, c0 : int
        Row and column of the first pixel of the ASIC
    arows, acols : int
        Number of rows and columns of the ASIC

    Returns
    -------
    nfs, nss : np.ndarray
        Unit vectors along the fast and slow scan directions
    """
    vfs = np.array((x[r0,c0+acols-1] - x[r0,c0],
                    y[r0,c0+acols-1] - y[r0,c0],
                    z[r0,c0+acols-1] - z[r0,c0]))
    vss = np.array((x[r0+arows-1,c0] - x[r0,c0],
                    y[r0+arows-1,c0] - y[r0,c0],
                    z[r0+arows-1,c0] - z[r0,c0]))
    return vfs/np.linalg.norm(vfs), vss/np.linalg.norm(vss)

"""
Dxtbx utility functions for writing detector models to experiment files and manipulating local frames.
"""


def write_experiment(detector, out_file):
    """
    Write a dxtbx Detector to a geometry file as a one-experiment ExperimentList

    dxtbx and DIALS read .expt files as an ExperimentList; a bare detector dictionary is not
    loadable, so the detector is wrapped in a single Experiment before being written.

    Parameters
    ----------
    detector : dxtbx.model.Detector
        Detector to serialize
    out_file : str
        Path to the output dxtbx .expt/.json geometry file
    """
    from dxtbx.model import Experiment, ExperimentList

    experiments = ExperimentList([Experiment(detector=detector)])
    experiments.as_file(out_file)
    print(f"Wrote dxtbx geometry file {out_file}.", flush=True)


def panel_layout(seg_geo):
    """
    Shape and first-pixel indices of the cctbx panels of one segment

    Jungfrau segments have double-width pixels at the inter-ASIC boundaries that span an
    integral number of pixel pitches, so cctbx represents the whole segment as one padded
    panel (512x1024 -> 514x1030). Other detectors, such as ePix, have fractional gaps that
    cannot be padded this way, and cctbx emits one panel per ASIC instead.

    Parameters
    ----------
    seg_geo : SegGeometry
        Segment geometry object from the segment geometry store

    Returns
    -------
    rows, cols : int, int
        Slow and fast size of each panel, in pixels
    indices : tuple
        One (row, col) pair per panel emitted for this segment
    """
    xx, yy = seg_geo.get_seg_xy_maps_um()
    pix = seg_geo.pixel_scale_size()
    # Extent of the segment in pixel pitches, counting its wide pixels
    cols = abs(xx[0, -1] - xx[0, 0]) / pix + 1
    rows = abs(yy[-1, 0] - yy[0, 0]) / pix + 1
    if abs(rows - round(rows)) < 1e-6 and abs(cols - round(cols)) < 1e-6:
        return int(round(rows)), int(round(cols)), ((0, 0),)
    return (*seg_geo.asic_rows_cols(), seg_geo.asic0indices())


def set_local_frame(node, basis):
    """
    Apply a cctbx basis to a dxtbx hierarchy node as its local frame

    Parameters
    ----------
    node : dxtbx.model.Panel or dxtbx.model.DetectorNode
        Node whose local frame is set
    basis : serialtbx.detector.basis
        Basis holding the orientation and translation of the node
    """
    from scitbx.matrix import col

    origin = col((basis * col((0, 0, 0, 1)))[0:3])
    fast = col((basis * col((1, 0, 0, 1)))[0:3]) - origin
    slow = col((basis * col((0, 1, 0, 1)))[0:3]) - origin
    node.set_local_frame(fast.elems, slow.elems, origin.elems)
"""
Tests for the pyFAI converter.

These cover the pixel corner array, which must straddle the psana pixel centers rather than
starting from them.
"""

import os

import numpy as np
import pytest

try:
    from LCLSGeom.converter import GeometryAccess, PsanaToPyFAI
except ImportError:
    pytest.skip("psana/PSCalib unavailable", allow_module_level=True)

from LCLSGeom.utils.frame import image_to_pyfai

TEMPLATES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src", "LCLSGeom", "templates")

DETECTORS = [
    ("jungfrau4M", 75.0),
    ("epix10k2M", 100.0),
]

# pyFAI stores the corner array as float32, so ~1e-7 relative precision on ~100 mm coordinates
PRECISION_MM = 1e-4


def template(detname):
    """Path to the shipped geometry template for a detector."""
    return os.path.join(TEMPLATES, f"geometry-def-{detname}.data")


def pyfai_frame_centers(detname):
    """psana pixel centers in mm, transformed the way get_pixel_corners transforms them."""
    geo = GeometryAccess(path=template(detname), pbits=0, use_wide_pix_center=False)
    camera = geo.get_top_geo().get_list_of_children()[0]
    camera.x0 = camera.y0 = camera.z0 = 0
    camera.rot_x = camera.rot_y = 0
    camera.rot_z -= 90
    x, y, z = camera.get_pixel_coords(do_tilt=True)
    X, Y, Z = image_to_pyfai(x, y, z, True)
    return X * 1e3, Y * 1e3, Z * 1e3


@pytest.mark.parametrize("detname,pixel_um", DETECTORS)
def test_pixel_corners_straddle_the_centers(detname, pixel_um):
    """
    The corner array extends half a pixel beyond the outermost pixel centers on every side.

    psana reports pixel centers, so corner offsets must be -1/2 and +1/2 of a pixel. Using 0
    and 1 instead puts the first corner on the center and translates the whole grid outward by
    half a pixel, which shows up as a lopsided bounding box: no margin on one side and a full
    pixel on the other.
    """
    detector = PsanaToPyFAI.convert(template(detname), detname)
    corners = detector.get_pixel_corners()
    x, y, _ = pyfai_frame_centers(detname)
    half = pixel_um * 1e-3 / 2

    # A half-pixel step along a tilted panel edge projects to slightly more than half a pixel,
    # so allow a few percent here. The bias this guards against is a full half pixel.
    tolerance = 0.05 * half
    for centers, corner in ((x, corners[..., 1] * 1e3), (y, corners[..., 2] * 1e3)):
        assert corner.min() - centers.min() == pytest.approx(-half, abs=tolerance)
        assert corner.max() - centers.max() == pytest.approx(half, abs=tolerance)
        # The grid must be centered on the pixel centers, not shifted off them
        drift = (corner.min() + corner.max()) / 2 - (centers.min() + centers.max()) / 2
        assert drift == pytest.approx(0.0, abs=PRECISION_MM)


@pytest.mark.parametrize("detname,pixel_um", DETECTORS)
def test_pixel_center_is_mean_of_its_corners(detname, pixel_um):
    """
    pyFAI derives a pixel center by averaging its four corners, which must return the psana
    center for that pixel.
    """
    detector = PsanaToPyFAI.convert(template(detname), detname)
    corners = detector.get_pixel_corners()
    x, y, _ = pyfai_frame_centers(detname)
    for row, col in [(0, 0), (0, 1), (1, 0), (5, 7)]:
        assert corners[row, col, :, 1].mean() * 1e3 == pytest.approx(x[0, row, col], abs=PRECISION_MM)
        assert corners[row, col, :, 2].mean() * 1e3 == pytest.approx(y[0, row, col], abs=PRECISION_MM)


@pytest.mark.parametrize("detname,pixel_um", DETECTORS)
def test_pixel_corners_span_one_pixel(detname, pixel_um):
    """Each pixel's corners enclose exactly one pixel pitch along both axes."""
    detector = PsanaToPyFAI.convert(template(detname), detname)
    corners = detector.get_pixel_corners()
    pixel_mm = pixel_um * 1e-3
    # Corners are ordered counterclockwise from the pixel's first corner
    for row, col in [(0, 0), (3, 11)]:
        quad = corners[row, col] * 1e3
        fast = np.linalg.norm(quad[1] - quad[0])
        slow = np.linalg.norm(quad[3] - quad[0])
        assert fast == pytest.approx(pixel_mm, abs=PRECISION_MM)
        assert slow == pytest.approx(pixel_mm, abs=PRECISION_MM)

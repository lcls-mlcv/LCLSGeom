"""
Tests for the CrystFEL converters.

The load-bearing one is the round trip: a .data file converted to .geom and back must land on the
same segment positions. That catches a swapped basis vector, which a one-way conversion does not.
"""

import os
import re

import numpy as np
import pytest

try:
    from LCLSGeom.converter import CrystFELToPsana, GeometryAccess, PsanaToCrystFEL
except ImportError:
    pytest.skip("psana/PSCalib unavailable", allow_module_level=True)

TEMPLATES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src", "LCLSGeom", "templates")

# template stem, psana detname -- these are different namespaces: the 16M template stem is
# jungfrau16M but psana calls the detector plain "jungfrau"
DETECTORS = [
    ("jungfrau05M", "jungfrau05M"),
    ("jungfrau1M", "jungfrau1M"),
    ("jungfrau4M", "jungfrau4M"),
    ("jungfrau16M", "jungfrau"),
    ("epix10k2M", "epix10k2M"),
]

# psana writes segment positions as whole microns, so a round trip cannot be exact
ROUND_TRIP_UM = 2.0


def template(stem):
    return os.path.join(TEMPLATES, f"geometry-def-{stem}.data")


def segment_centers(path):
    """Mean x, y of every segment of a psana geometry file, in microns."""
    geo = GeometryAccess(path=path, pbits=0, use_wide_pix_center=False)
    seg = geo.get_seg_geo().algo
    x, y, _ = geo.get_pixel_coords(do_tilt=True)
    nsegs = int(x.size / seg.size())
    x = x.reshape(nsegs, *seg.shape())
    y = y.reshape(nsegs, *seg.shape())
    return np.array([[x[n].mean(), y[n].mean()] for n in range(nsegs)])


@pytest.mark.parametrize("stem,detname", DETECTORS)
def test_round_trip_preserves_segment_positions(tmp_path, stem, detname):
    """
    .data -> .geom -> .data lands back on the same segment positions.

    A basis vector used in place of its partner shifts segments by millimetres while leaving the
    file structurally valid, so only a position comparison catches it.
    """
    geom = str(tmp_path / f"{stem}.geom")
    data = str(tmp_path / f"{stem}.data")
    PsanaToCrystFEL.convert(in_file=template(stem), detname=detname, out_file=geom)
    CrystFELToPsana.convert(in_file=geom, detname=detname, out_file=data)

    before, after = segment_centers(template(stem)), segment_centers(data)
    assert before.shape == after.shape
    assert np.abs(before - after).max() < ROUND_TRIP_UM


@pytest.mark.parametrize("stem,detname", DETECTORS)
def test_group_definitions_appear_once(tmp_path, stem, detname):
    """
    The millepede groups are a per-file trailer, not a per-ASIC one.

    They used to be emitted inside the ASIC loop, repeating the whole block once per ASIC (64
    times for a Jungfrau4M) and tripling the file size.
    """
    geom = str(tmp_path / f"{stem}.geom")
    PsanaToCrystFEL.convert(in_file=template(stem), detname=detname, out_file=geom)
    lines = [l for l in open(geom).read().splitlines() if l.startswith("group_")]

    assert lines, "expected group definitions for a segmented detector"
    assert len(lines) == len(set(lines))


@pytest.mark.parametrize("stem,detname", DETECTORS)
def test_every_segment_is_written(tmp_path, stem, detname):
    """
    Each psana segment reaches the .geom and comes back in the .data.

    The 16M panel list used to stop at 16 of its 32 segments, so half the detector was silently
    dropped on the way back.
    """
    geom = str(tmp_path / f"{stem}.geom")
    data = str(tmp_path / f"{stem}.data")
    PsanaToCrystFEL.convert(in_file=template(stem), detname=detname, out_file=geom)
    CrystFELToPsana.convert(in_file=geom, detname=detname, out_file=data)

    nsegs = len(segment_centers(template(stem)))
    panels = {int(p) for p in re.findall(r"^p(\d+)a\d+/fs", open(geom).read(), re.M)}
    assert panels == set(range(nsegs))

    # One record per segment, plus the IP-detector record that carries the distance
    records = [l for l in open(data).read().splitlines() if l.strip() and not l.startswith("#")]
    assert len(records) == nsegs + 1

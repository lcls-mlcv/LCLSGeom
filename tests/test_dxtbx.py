"""
Tests for the dxtbx converters, PsanaToDXTBX and PsanaToCCTBX.

Both converters are checked against the raw psana pixel coordinates of the shipped geometry
templates, so no calibration database access is required. That comparison is the load-bearing
one: a unit error shows up as a factor of 1000, a frame error as a reflection or axis swap, and
a center/corner error as a clean half-pixel offset.
"""

import json
import os

import numpy as np
import pytest

pytest.importorskip("dxtbx")
try:
    from LCLSGeom.converter import PsanaToCCTBX, PsanaToDXTBX
except ImportError:
    pytest.skip("psana/PSCalib unavailable", allow_module_level=True)

from dxtbx.model import DetectorFactory
from dxtbx.model.experiment_list import ExperimentListFactory

from LCLSGeom.utils.pscalib_utils import get_root_geometry

TEMPLATES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src", "LCLSGeom", "templates")

# detname, calib shape, pixel size in microns, ASIC layout
DETECTORS = [
    ("jungfrau4M", (8, 512, 1024), 75.0, (2, 4)),
    ("epix10k2M", (16, 352, 384), 100.0, (2, 2)),
]

# Sub-pixel tolerance in mm: well under one pixel for either detector
TOLERANCE_MM = 1e-3


def template(detname):
    """Path to the shipped geometry template for a detector."""
    return os.path.join(TEMPLATES, f"geometry-def-{detname}.data")


def psana_coords(detname, shape, zeroed):
    """
    Pixel center coordinates in mm, as psana reports them.

    ``zeroed`` is passed straight through as ``image_frame``, so the coordinates come out in the
    same frame the converter built its panels in and are directly comparable to them.
    """
    converter = PsanaToDXTBX(template(detname), detname)
    camera = get_root_geometry(geo=converter.geo, image_frame=zeroed)
    x, y, z = camera.get_pixel_coords(do_tilt=True)
    x.shape = y.shape = z.shape = shape
    return x * 1e-3, y * 1e-3, z * 1e-3


def probe_pixels(arows, acols):
    """Corner and center pixels of an ASIC, as (slow, fast) offsets."""
    return [(0, 0), (0, acols - 1), (arows - 1, 0), (arows - 1, acols - 1), (arows // 2, acols // 2)]


@pytest.mark.parametrize("detname,shape", [(d[0], d[1]) for d in DETECTORS])
def test_dxtbx_pixel_positions_match_psana(detname, shape):
    """Every probed pixel center agrees with psana to well under a pixel."""
    detector = PsanaToDXTBX.convert(template(detname), detname)
    x, y, z = psana_coords(detname, shape, zeroed=True)
    converter = PsanaToDXTBX(template(detname), detname)
    seg = converter.geo.get_seg_geo().algo
    arows, acols = seg.asic_rows_cols()
    indices = seg.asic0indices()

    errors = []
    for i, panel in enumerate(detector):
        n, a = divmod(i, len(indices))
        r0, c0 = indices[a]
        for dr, dc in probe_pixels(arows, acols):
            # dxtbx indexes pixel centers at half-integer coordinates
            lab = np.array(panel.get_pixel_lab_coord((dc + 0.5, dr + 0.5)))
            ref = np.array((x[n, r0 + dr, c0 + dc], y[n, r0 + dr, c0 + dc], z[n, r0 + dr, c0 + dc]))
            errors.append(np.linalg.norm(lab - ref))
    assert max(errors) < TOLERANCE_MM


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_panels_are_per_asic_and_placement_is_zeroed(detname, shape, pixel_um, asics):
    """
    One panel per ASIC, each the size of an ASIC, with the detector placement removed.

    Zeroing removes the placement, not the per-module metrology, so a detector whose modules
    are tilted keeps small out-of-plane offsets. What must go is the sample-detector distance:
    the panels straddle z = 0 rather than sitting at -detz.
    """
    detector = PsanaToDXTBX.convert(template(detname), detname)
    nsegments, rows, cols = shape
    assert len(detector) == nsegments * asics[0] * asics[1]
    for panel in detector:
        assert panel.get_image_size() == (cols // asics[1], rows // asics[0])
        assert panel.get_pixel_size() == pytest.approx((pixel_um * 1e-3, pixel_um * 1e-3))

    z = np.array([panel.get_origin()[2] for panel in detector])
    assert np.mean(z) == pytest.approx(0.0, abs=0.01)
    assert np.max(np.abs(z)) < 1.0


@pytest.mark.parametrize("detname", [d[0] for d in DETECTORS])
def test_dxtbx_serialized_origins_are_absolute(detname):
    """
    Serialized panel origins are absolute, not relative to a parent group.

    dxtbx always writes a hierarchy, but for this converter its root is the identity and every
    panel hangs directly off it. Consumers that snapshot ``to_dict()["panels"]`` and rebuild
    without the hierarchy therefore keep the correct positions; a nested model would collapse.
    """
    detector = PsanaToDXTBX.convert(template(detname), detname)
    as_dict = detector.to_dict()
    root = as_dict["hierarchy"]
    assert len(root["children"]) == len(detector)
    assert root["origin"] == pytest.approx((0.0, 0.0, 0.0))

    rebuilt = DetectorFactory.from_dict({"panels": as_dict["panels"]})
    for before, after in zip(detector, rebuilt):
        assert np.allclose(before.get_origin(), after.get_origin())


@pytest.mark.parametrize("detname", [d[0] for d in DETECTORS])
def test_detname_must_match_the_file(detname):
    """A detector name that disagrees with the file's segment count is a caller error."""
    with pytest.raises(ValueError):
        PsanaToDXTBX.convert(template(detname), "rayonix")


def test_cctbx_pixel_positions_match_psana():
    """
    Padded panel pixels agree with psana once the inserted gap pixels are accounted for.

    Mapping a psana (row, col) to its padded index reproduces the psana coordinate, including
    at the double-width boundaries. cctbx uses the first pixel center as the panel origin, so
    pixel centers sit at integer coordinates here rather than half-integer ones.
    """
    detname, shape = "jungfrau4M", (8, 512, 1024)
    detector = PsanaToCCTBX(template(detname), detname).build_detector()
    x, y, z = psana_coords(detname, shape, zeroed=False)

    def padded_index(row, col):
        return row + (2 if row >= 256 else 0), col + 2 * sum(1 for b in (256, 512, 768) if col >= b)

    errors = []
    for n, panel in enumerate(detector):
        for r, c in [(0, 0), (0, 1023), (511, 0), (511, 1023), (255, 255), (256, 256), (300, 700)]:
            slow, fast = padded_index(r, c)
            lab = np.array(panel.get_pixel_lab_coord((fast, slow)))
            ref = np.array((x[n, r, c], y[n, r, c], z[n, r, c]))
            errors.append(np.linalg.norm(lab[:2] - ref[:2]))
    assert max(errors) < TOLERANCE_MM


def test_cctbx_jungfrau_panels_are_padded():
    """
    A Jungfrau segment becomes one panel padded for its double-width pixels, in a hierarchy.

    The inter-ASIC pixels span three pixel pitches, so a 512x1024 segment spans 514x1030
    pitches and cctbx represents it as a single panel of that size. Unlike PsanaToDXTBX,
    nothing is zeroed and the panel origins are local to their module group.
    """
    detector = PsanaToCCTBX(template("jungfrau4M"), "jungfrau4M").build_detector()
    assert len(detector) == 8
    for panel in detector:
        assert panel.get_image_size() == (1030, 514)
    assert detector[0].get_distance() > 0

    as_dict = detector.to_dict()
    assert set(as_dict.keys()) == {"panels", "hierarchy"}
    root = as_dict["hierarchy"]
    assert root["origin"][2] < 0
    assert len(root["children"]) == 8
    for group in root["children"]:
        assert [child["panel"] for child in group["children"]]

    local = as_dict["panels"][0]["origin"]
    assert local[0] == pytest.approx(-(1030 - 1) / 2 * 0.075)
    assert abs(local[1]) == pytest.approx((514 - 1) / 2 * 0.075)
    assert local[2] == pytest.approx(0.0)


def test_cctbx_epix_panels_are_per_asic():
    """
    An ePix segment cannot be padded, so cctbx emits one panel per ASIC.

    Its inter-ASIC gaps span a fractional number of pixel pitches, unlike Jungfrau's, so the
    padded-versus-per-ASIC split has to be derived from the segment map rather than hardcoded.
    """
    detector = PsanaToCCTBX(template("epix10k2M"), "epix10k2M").build_detector()
    assert len(detector) == 16 * 4
    for panel in detector:
        assert panel.get_image_size() == (192, 176)


@pytest.mark.parametrize("write,panels", [
    (lambda i, d, o: PsanaToDXTBX.convert_to_file(i, d, o), 64),
    (lambda i, d, o: PsanaToCCTBX.convert_to_file(i, d, o), 8),
])
def test_written_geometry_is_loadable_by_dials(tmp_path, write, panels):
    """
    The written file is a real .expt, i.e. an ExperimentList that DIALS can open.

    A bare detector dictionary is not loadable, so the detector has to be wrapped in an
    Experiment; this fails with InvalidExperimentListError otherwise.
    """
    out_file = tmp_path / "geometry.expt"
    write(template("jungfrau4M"), "jungfrau4M", str(out_file))

    written = json.loads(out_file.read_text())
    assert written["__id__"] == "ExperimentList"
    assert len(written["experiment"]) == 1

    experiments = ExperimentListFactory.from_json_file(str(out_file), check_format=False)
    assert len(experiments) == 1
    assert len(experiments[0].detector) == panels


def test_real_geometry_from_calibration(tmp_path):
    """
    Smoke test over a real metrology file, which needs a psana environment.

    The MFX detector is a Jungfrau 16M, whose psana name is plain "jungfrau": 32 segments, so
    256 ASIC panels or 32 padded module panels.
    """
    from LCLSGeom.manager import get_geometry

    try:
        geometry_file = get_geometry("jungfrau", "mfx101591226", 20)
    except (Exception, SystemExit) as exc:
        # psana exits rather than raising when the experiment's xtc files are unreachable
        pytest.skip(f"calibration lookup unavailable: {exc}")

    assert len(PsanaToDXTBX.convert(geometry_file, "jungfrau")) == 256
    out_file = tmp_path / "cctbx.expt"
    PsanaToCCTBX.convert_to_file(geometry_file, "jungfrau", str(out_file))
    experiments = ExperimentListFactory.from_json_file(str(out_file), check_format=False)
    assert len(experiments[0].detector) == 32

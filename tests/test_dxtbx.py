"""
Tests for the dxtbx converters, PsanaToDXTBX and PsanaToCCTBX.

Both converters are checked against the raw psana pixel coordinates of the shipped geometry
templates, so no calibration database access is required.
"""

import json
import os

import numpy as np
import pytest

pytest.importorskip("dxtbx")
try:
    from LCLSGeom.converter import PsanaToCCTBX, PsanaToDXTBX, PsanaToPyFAI
except ImportError:
    pytest.skip("psana/PSCalib unavailable", allow_module_level=True)

from dxtbx.model import DetectorFactory

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

    With ``zeroed`` the camera placement is removed the same way PsanaToDXTBX removes it, so
    the result is directly comparable to that converter's panels.
    """
    converter = PsanaToDXTBX(template(detname), detname)
    if zeroed:
        camera = converter.zero_placement()
    else:
        camera = converter.geo.get_top_geo().get_list_of_children()[0]
    x, y, z = camera.get_pixel_coords(do_tilt=True)
    x.shape = y.shape = z.shape = shape
    return x * 1e-3, y * 1e-3, z * 1e-3


def probe_pixels(arows, acols):
    """Corner and center pixels of an ASIC, as (slow, fast) offsets."""
    return [(0, 0), (0, acols - 1), (arows - 1, 0), (arows - 1, acols - 1), (arows // 2, acols // 2)]


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_panel_count_and_shape(detname, shape, pixel_um, asics):
    """One panel per ASIC, each the size of an ASIC, with the expected pixel pitch."""
    detector = PsanaToDXTBX.convert(template(detname), detname)
    nsegments, rows, cols = shape
    nasics = asics[0] * asics[1]
    assert len(detector) == nsegments * nasics
    for panel in detector:
        assert panel.get_image_size() == (cols // asics[1], rows // asics[0])
        assert panel.get_pixel_size() == pytest.approx((pixel_um * 1e-3, pixel_um * 1e-3))


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_pixel_positions_match_psana(detname, shape, pixel_um, asics):
    """
    Every probed pixel center agrees with psana to well under a pixel.

    This is the load-bearing test: a unit error shows up as a factor of 1000, a frame error as
    a reflection or axis swap, and a center/corner error as a clean half-pixel offset.
    """
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
def test_dxtbx_axes_are_orthonormal(detname, shape, pixel_um, asics):
    """Panel basis vectors are unit length and mutually perpendicular."""
    detector = PsanaToDXTBX.convert(template(detname), detname)
    for panel in detector:
        fast = np.array(panel.get_fast_axis())
        slow = np.array(panel.get_slow_axis())
        assert np.linalg.norm(fast) == pytest.approx(1.0)
        assert np.linalg.norm(slow) == pytest.approx(1.0)
        assert np.dot(fast, slow) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_placement_is_zeroed(detname, shape, pixel_um, asics):
    """
    Panel origins hold metrology only, with the distance available separately.

    A detector at zero distance is degenerate for dxtbx, so the beam center cannot be
    computed. That is expected here and is pinned so it cannot change silently.
    """
    converter = PsanaToDXTBX(template(detname), detname)
    converter.get_panels()
    assert converter.detz_mm > 0

    # Zeroing removes the detector placement, not the per-module metrology, so the panels of a
    # detector whose modules are tilted keep their small out-of-plane offsets. What must go is
    # the sample-detector distance: the panels straddle z = 0 instead of sitting at -detz.
    detector = PsanaToDXTBX.convert(template(detname), detname)
    z = np.array([panel.get_origin()[2] for panel in detector])
    assert np.mean(z) == pytest.approx(0.0, abs=0.01)
    assert np.max(np.abs(z)) < 1.0


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_serialized_origins_are_absolute(detname, shape, pixel_um, asics):
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
    for serialized, panel in zip(as_dict["panels"], detector):
        assert np.allclose(serialized["origin"], panel.get_origin())

    # Rebuilding from panels alone, the way BraggFAI does, preserves the geometry
    rebuilt = DetectorFactory.from_dict({"panels": as_dict["panels"]})
    for before, after in zip(detector, rebuilt):
        assert np.allclose(before.get_origin(), after.get_origin())


def test_dxtbx_planar_detector_has_no_beam_centre():
    """
    A planar detector ends up exactly at zero distance, which dxtbx cannot intersect.

    This is the documented consequence of zeroing the placement, pinned here so it cannot
    change silently. Callers needing an absolutely placed model should use PsanaToCCTBX.
    """
    detector = PsanaToDXTBX.convert(template("jungfrau4M"), "jungfrau4M")
    assert all(panel.get_origin()[2] == pytest.approx(0.0) for panel in detector)
    with pytest.raises(RuntimeError):
        detector[0].get_beam_centre((0.0, 0.0, -1.0))


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_json_round_trip(detname, shape, pixel_um, asics):
    """
    Serializing and reloading preserves the geometry.

    dxtbx re-orthonormalizes the basis vectors on load, so values are compared numerically
    rather than for exact equality.
    """
    detector = PsanaToDXTBX.convert(template(detname), detname)
    reloaded = DetectorFactory.from_dict(json.loads(json.dumps(detector.to_dict())))
    assert len(reloaded) == len(detector)
    for before, after in zip(detector, reloaded):
        assert np.allclose(before.get_origin(), after.get_origin())
        assert np.allclose(before.get_fast_axis(), after.get_fast_axis())
        assert np.allclose(before.get_slow_axis(), after.get_slow_axis())
        assert before.get_image_size() == after.get_image_size()


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_inter_panel_distances_are_rigid(detname, shape, pixel_um, asics):
    """
    Distances between panel origins match the same distances computed from psana.

    A single-panel test cannot catch a per-panel placement error; this can.
    """
    detector = PsanaToDXTBX.convert(template(detname), detname)
    x, y, z = psana_coords(detname, shape, zeroed=True)
    converter = PsanaToDXTBX(template(detname), detname)
    seg = converter.geo.get_seg_geo().algo
    indices = seg.asic0indices()

    # Compare pixel centers on both sides: panel origins are outer corners, and that
    # half-pixel shift does not cancel between panels that are tilted differently.
    origins = np.array([panel.get_pixel_lab_coord((0.5, 0.5)) for panel in detector])
    reference = np.array([
        (x[n, r0, c0], y[n, r0, c0], z[n, r0, c0])
        for n in range(shape[0])
        for (r0, c0) in indices
    ])
    for i in range(0, len(origins), 7):
        for j in range(i + 1, len(origins), 11):
            got = np.linalg.norm(origins[i] - origins[j])
            want = np.linalg.norm(reference[i] - reference[j])
            assert got == pytest.approx(want, abs=TOLERANCE_MM)


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_agrees_with_pyfai(detname, shape, pixel_um, asics):
    """
    Cross-check pixel centers against the pyFAI converter.

    Both converters describe the same detector, so their pixel centers must coincide. This
    catches a unit error as a factor of 1000 and a frame error as a reflection or axis swap.
    """
    detector = PsanaToDXTBX.convert(template(detname), detname)
    pyfai = PsanaToPyFAI.convert(template(detname), detname)
    xc, yc, zc = pyfai.calc_cartesian_positions()
    converter = PsanaToDXTBX(template(detname), detname)
    seg = converter.geo.get_seg_geo().algo
    indices = seg.asic0indices()
    ss, fs = pyfai.ss_size, pyfai.fs_size
    nasics_slow, nasics_fast = pyfai.asics_shape

    errors = []
    for i, panel in enumerate(detector):
        n, a = divmod(i, len(indices))
        arow, acol = a // nasics_fast, a % nasics_fast
        for dr, dc in probe_pixels(ss, fs):
            row = n * nasics_slow * ss + arow * ss + dr
            col = acol * fs + dc
            # pyFAI (x, y, z) is the psana image frame as (y, x, -z); see utils.frame
            ref = np.array([yc[row, col], xc[row, col], -zc[row, col]]) * 1e3
            lab = np.array(panel.get_pixel_lab_coord((dc + 0.5, dr + 0.5)))
            errors.append(np.linalg.norm(ref[:2] - lab[:2]))
    assert max(errors) < TOLERANCE_MM


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_dxtbx_detname_must_match_the_file(detname, shape, pixel_um, asics):
    """A detector name that disagrees with the file's segment count is a caller error."""
    with pytest.raises(ValueError):
        PsanaToDXTBX.convert(template(detname), "rayonix")


@pytest.mark.parametrize("converter,panels", [(PsanaToDXTBX, 64), (PsanaToCCTBX, 8)])
def test_convert_to_file_is_loadable_by_dials(tmp_path, converter, panels):
    """
    The written file is a real .expt, i.e. an ExperimentList that DIALS can open.

    A bare detector dictionary is not loadable, so the detector has to be wrapped in an
    Experiment; this fails with InvalidExperimentListError otherwise.
    """
    from dxtbx.model.experiment_list import ExperimentListFactory

    out_file = tmp_path / "geometry.expt"
    converter.convert_to_file(template("jungfrau4M"), "jungfrau4M", str(out_file))

    written = json.loads(out_file.read_text())
    assert written["__id__"] == "ExperimentList"
    assert len(written["experiment"]) == 1

    experiments = ExperimentListFactory.from_json_file(str(out_file), check_format=False)
    assert len(experiments) == 1
    assert len(experiments[0].detector) == panels


@pytest.mark.parametrize("detname,shape,pixel_um,asics", DETECTORS)
def test_cctbx_hierarchy_structure(detname, shape, pixel_um, asics):
    """
    The cctbx model is a three level hierarchy with one group per psana segment.

    Panel origins are local to their group, and the sample-detector distance is kept rather
    than zeroed.
    """
    detector = PsanaToCCTBX.convert(template(detname), detname)
    as_dict = detector.to_dict()
    assert set(as_dict.keys()) == {"panels", "hierarchy"}

    root = as_dict["hierarchy"]
    assert len(root["children"]) == shape[0]
    assert root["origin"][2] < 0
    for group in root["children"]:
        for child in group["children"]:
            assert "panel" in child

    assert detector[0].get_distance() > 0
    for panel in detector:
        assert panel.get_pixel_size() == pytest.approx((pixel_um * 1e-3, pixel_um * 1e-3))


def test_cctbx_jungfrau_panels_are_padded():
    """
    A Jungfrau segment becomes one panel padded for its double-width pixels.

    The inter-ASIC pixels span three pixel pitches, so a 512x1024 segment spans 514x1030
    pitches and cctbx represents it as a single panel of that size.
    """
    detector = PsanaToCCTBX.convert(template("jungfrau4M"), "jungfrau4M")
    assert len(detector) == 8
    for panel in detector:
        assert panel.get_image_size() == (1030, 514)

    local = detector.to_dict()["panels"][0]["origin"]
    assert local[0] == pytest.approx(-(1030 - 1) / 2 * 0.075)
    assert abs(local[1]) == pytest.approx((514 - 1) / 2 * 0.075)
    assert local[2] == pytest.approx(0.0)


def test_cctbx_epix_panels_are_per_asic():
    """
    An ePix segment cannot be padded, so cctbx emits one panel per ASIC.

    Its inter-ASIC gaps span a fractional number of pixel pitches, unlike Jungfrau's.
    """
    detector = PsanaToCCTBX.convert(template("epix10k2M"), "epix10k2M")
    assert len(detector) == 16 * 4
    for panel in detector:
        assert panel.get_image_size() == (192, 176)


def test_cctbx_pixel_positions_match_psana():
    """
    Padded panel pixels agree with psana once the inserted gap pixels are accounted for.

    Mapping a psana (row, col) to its padded index reproduces the psana coordinate exactly,
    including at the double-width boundaries. Note that cctbx uses the first pixel center as
    the panel origin, so pixel centers sit at integer coordinates here.
    """
    detname, shape = "jungfrau4M", (8, 512, 1024)
    detector = PsanaToCCTBX.convert(template(detname), detname)
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


def test_cctbx_json_round_trip():
    """Reloading a serialized cctbx model preserves absolute panel positions."""
    detector = PsanaToCCTBX.convert(template("jungfrau4M"), "jungfrau4M")
    reloaded = DetectorFactory.from_dict(json.loads(json.dumps(detector.to_dict())))
    assert len(reloaded) == len(detector)
    for before, after in zip(detector, reloaded):
        assert np.allclose(before.get_origin(), after.get_origin())


def test_cctbx_matches_reference_expt():
    """
    The emitted structure matches a refined geometry produced by cctbx at MFX.

    Origins differ, since the reference is refined and ours comes from a template, so only the
    structure and panel metadata are compared.
    """
    reference_file = "/sdf/data/lcls/ds/mfx/mfx101817926/results/common/geom/r25_bayfaigeom.expt"
    if not os.access(reference_file, os.R_OK):
        pytest.skip("reference geometry not readable")
    with open(reference_file) as f:
        reference = json.load(f)["detector"][0]

    # "jungfrau" is the psana name for the 16M, which is what MFX ran here
    detector = PsanaToCCTBX.convert(template("jungfrau16M"), "jungfrau").to_dict()
    assert len(detector["panels"]) == len(reference["panels"])
    assert len(detector["hierarchy"]["children"]) == len(reference["hierarchy"]["children"])
    for ours, theirs in zip(detector["panels"], reference["panels"]):
        assert tuple(ours["image_size"]) == tuple(theirs["image_size"])
        assert ours["pixel_size"] == pytest.approx(theirs["pixel_size"])
        # Panel origins are local to their module group in both models
        assert np.allclose(np.abs(ours["origin"]), np.abs(theirs["origin"]))


@pytest.mark.parametrize("converter,panels", [(PsanaToDXTBX, 256), (PsanaToCCTBX, 32)])
def test_real_geometry_from_calibration(converter, panels):
    """
    Smoke test over a real metrology file, which needs a psana environment.

    The MFX detector is a Jungfrau 16M, whose psana name is plain "jungfrau": 32 segments, so
    256 ASIC panels or 32 padded module panels.
    """
    from LCLSGeom.manager import get_geometry

    try:
        geometry_file = get_geometry("jungfrau", "mfx101591226", 20)
    except Exception as exc:
        pytest.skip(f"calibration lookup unavailable: {exc}")
    detector = converter.convert(geometry_file, "jungfrau")
    assert len(detector) == panels

"""
Tests for the detector name conventions and the generated CrystFEL group definitions.

These need no psana environment: the tables and generators are plain data, which is the point of
keeping them free of heavy imports.
"""

import pytest

from LCLSGeom.utils.crystfel_utils import panel_groups, panel_names, tailer_crystfel
from LCLSGeom.utils.detectors import FAMILIES, get_family, get_layout
from LCLSGeom.utils.pscalib_utils import (PSANA_DETECTORS, get_calib_detname, get_full_detname,
                                          group_from_detname, source_from_detname)

# psana detname, family, template stem, calibration source
DETECTORS = [
    ("Rayonix",       "Rayonix",      "rayonix",       "Rayonix"),
    ("jungfrau05M",   "Jungfrau05M",  "jungfrau05M",   "Jungfrau"),
    ("jungfrau1M",    "Jungfrau1M",   "jungfrau1M",    "Jungfrau"),
    ("jungfrau4M",    "Jungfrau4M",   "jungfrau4M",    "Jungfrau"),
    # psana calls the 16M plain "jungfrau"; its template stem does not
    ("jungfrau",      "Jungfrau16M",  "jungfrau16M",   "Jungfrau"),
    ("epix10k2M",     "ePix10k2M",    "epix10k2M",     "Epix10ka2M"),
    ("Epix10kaQuad0", "ePix10kaQuad", "epix10kaQuad0", "Epix10kaQuad.0"),
    ("Epix10kaQuad3", "ePix10kaQuad", "epix10kaQuad3", "Epix10kaQuad.3"),
]


@pytest.mark.parametrize("detname,family,stem,source", DETECTORS)
def test_name_conventions(detname, family, stem, source):
    """Each software's spelling of a detector is reachable from the psana name."""
    assert get_family(detname) == family
    assert get_full_detname(detname) == stem
    assert get_calib_detname(detname) == source


@pytest.mark.parametrize("detname", [d[0] for d in DETECTORS])
def test_unknown_names_raise(detname):
    """A name that is not a psana detector is a caller error, in every mapping."""
    for fn in (get_family, get_full_detname, get_calib_detname):
        fn(detname)  # the real name works
        with pytest.raises(ValueError):
            fn(detname + "-nonesuch")


def test_source_carries_the_quad_index_once():
    """
    The ePix10ka quad index belongs to the source but not to the calibration group.

    A quad source must not pick up a trailing '.0' on top of its own index, and all four quads
    share one group.
    """
    assert source_from_detname("Epix10kaQuad2", "cxi") == "CxiDs1.0:Epix10kaQuad.2"
    assert source_from_detname("jungfrau4M", "cxi") == "CxiDs1.0:Jungfrau.0"
    assert group_from_detname("Epix10kaQuad2") == "Epix10kaQuad::CalibV1"
    assert group_from_detname("Epix10kaQuad0") == group_from_detname("Epix10kaQuad3")

    with pytest.raises(ValueError):
        source_from_detname("jungfrau4M", "nonesuch")


@pytest.mark.parametrize("detname", list(PSANA_DETECTORS))
def test_panel_names_cover_every_segment(detname):
    """One CrystFEL panel per psana segment, in segment order."""
    names = panel_names(detname)
    if names is None:
        assert get_family(detname) == "Rayonix"
        return
    nsegs, _ = get_layout(get_family(detname))
    assert names.split(",") == [f"p{n}a0" for n in range(nsegs)]


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_panel_groups_partition_the_asics(family):
    """
    Every ASIC appears in exactly one panel group, and every segment in at most one quadrant.

    This is what millepede needs: the groups are a partition, not an overlapping selection.
    """
    groups = panel_groups(family)
    if not groups:
        assert family == "Rayonix"
        return
    nsegs, nasics = get_layout(family)

    asics = [a for n in range(nsegs) for a in groups[f"group_p{n}"].split(",")]
    assert sorted(asics) == sorted(f"p{n}a{a}" for n in range(nsegs) for a in range(nasics))

    quads = [k for k in groups if k.startswith("group_q")]
    if quads:
        assert len(quads) == 4
        assert groups["group_all"] == "q0,q1,q2,q3"
        members = [p for q in quads for p in groups[q].split(",")]
        assert sorted(members, key=lambda s: int(s[1:])) == [f"p{n}" for n in range(nsegs)]
    else:
        # Fewer than two segments per quadrant: a quadrant would not be a useful rigid body
        assert nsegs % 4 != 0 or nsegs // 4 < 2


def test_tailer_is_empty_for_a_monolithic_detector():
    """Rayonix has no panels to group, so its tailer carries no group definitions."""
    assert tailer_crystfel("Rayonix").strip() == ""
    assert tailer_crystfel("Jungfrau16M").count("group_p") == 32

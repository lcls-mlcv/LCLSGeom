"""
Physical facts about the LCLS detectors, and the shapes each one is allowed to have.
"""

psana_detnames = (
    "UNDEFINED",
    "Rayonix",
    "jungfrau05M",
    "jungfrau1M",
    "jungfrau4M",
    "jungfrau",
    "epix10k2M",
    "Epix10kaQuad0",
    "Epix10kaQuad1",
    "Epix10kaQuad2",
    "Epix10kaQuad3",
)

# The detector family, one per psana name above. Several psana names collapse onto one family:
# the four ePix10ka quads are all 'ePix10kaQuad', and psana's plain 'jungfrau' is the 16M.
family_names = (
    "UNDEFINED",
    "Rayonix",
    "Jungfrau05M",
    "Jungfrau1M",
    "Jungfrau4M",
    "Jungfrau16M",
    "ePix10k2M",
    "ePix10kaQuad",
    "ePix10kaQuad",
    "ePix10kaQuad",
    "ePix10kaQuad",
)

psana_detnames_lower = tuple(name.lower() for name in psana_detnames)

psana_to_family = dict(zip(psana_detnames_lower, family_names))


def get_family(detname: str) -> str:
    """
    Start from a psana detector name and return its detector family.

    Parameters
    ----------
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)

    Returns
    -------
    family : str
        Detector family (e.g., 'Jungfrau4M', 'ePix10k2M', etc.)
    """
    family = psana_to_family.get(detname.lower(), "UNDEFINED")
    if family == "UNDEFINED":
        raise ValueError(f"Unknown detector type: {detname}")
    return family


# Family, pixel size [um], calib shape, ASIC layout (slow, fast) per segment.
DETECTOR_SPECS = {
    "Rayonix":      ("Rayonix",       44.0, (1, 7680, 7680), (1, 1)),
    "Rayonix2x2":   ("Rayonix",       88.0, (1, 3840, 3840), (1, 1)),
    "Rayonix3x3":   ("Rayonix",      132.0, (1, 2560, 2560), (1, 1)),
    "Rayonix4x4":   ("Rayonix",      176.0, (1, 1920, 1920), (1, 1)),
    "Rayonix5x5":   ("Rayonix",      220.0, (1, 1536, 1536), (1, 1)),
    "Rayonix6x6":   ("Rayonix",      264.0, (1, 1280, 1280), (1, 1)),
    "Rayonix8x8":   ("Rayonix",      352.0, (1,  960,  960), (1, 1)),
    "Rayonix10x10": ("Rayonix",      440.0, (1,  768,  768), (1, 1)),
    "ePix10k2M":    ("ePix10k2M",    100.0, (16, 352, 384),  (2, 2)),
    "ePix10kaQuad": ("ePix10kaQuad", 100.0, (4,  352, 384),  (2, 2)),
    "Jungfrau05M":  ("Jungfrau05M",   75.0, (1,  512, 1024), (2, 4)),
    "Jungfrau1M":   ("Jungfrau1M",    75.0, (2,  512, 1024), (2, 4)),
    "Jungfrau4M":   ("Jungfrau4M",    75.0, (8,  512, 1024), (2, 4)),
    "Jungfrau16M":  ("Jungfrau16M",   75.0, (32, 512, 1024), (2, 4)),
}

FAMILIES = frozenset(spec[0] for spec in DETECTOR_SPECS.values())

SUPPORTED_SHAPES = {
    (family, shape) for family, _, shape, _ in DETECTOR_SPECS.values()
}


def get_layout(family: str) -> tuple:
    """
    Number of segments and of ASICs per segment, for a detector family.

    Parameters
    ----------
    family : str
        Detector family, as returned by get_family

    Returns
    -------
    nsegs, nasics : int, int
        Segments in the detector, and ASICs in one segment
    """
    for spec_family, _, shape, asics in DETECTOR_SPECS.values():
        if spec_family == family:
            return shape[0], asics[0] * asics[1]
    raise ValueError(f"Unknown detector family: {family}")


def check_detname(detname: str, shape: tuple) -> None:
    """
    Raise if a detector name is inconsistent with the shape of its geometry file.

    Keyed on the detector family, so a dxtbx, cctbx or CrystFEL caller can validate a shape
    without routing through any one consumer's naming convention.

    Parameters
    ----------
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
    shape : tuple
        Shape of the geometry file as (<number-of-segments>, <rows>, <cols>)
    """
    key = (get_family(detname), shape)
    if key not in SUPPORTED_SHAPES:
        raise ValueError(
            f"Detector not recognized: detname={detname!r}, shape={shape}. "
            f"Supported (family, shape) pairs: {sorted(SUPPORTED_SHAPES)}"
        )

from pyFAI.detectors import Detector
from LCLSGeom.utils.calibmap import get_pyFAI_detname

class Rayonix(Detector):
    """
    PyFAI Detector instance for the Rayonix
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000044
        self.pixel_size_um = 44.0
        self.detname = "Rayonix"
        self.calib_shape = (1, 7680, 7680)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 1
        self.asics_shape = (1, 1)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Rayonix2x2(Detector):
    """
    PyFAI Detector instance for the Rayonix 2x2 binning
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000088
        self.pixel_size_um = 88.0
        self.detname = "Rayonix"
        self.calib_shape = (1, 3840, 3840)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 1
        self.asics_shape = (1, 1)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Rayonix3x3(Detector):
    """
    PyFAI Detector instance for the Rayonix 3x3 binning
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000132
        self.pixel_size_um = 132.0
        self.detname = "Rayonix"
        self.calib_shape = (1, 2560, 2560)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 1
        self.asics_shape = (1, 1)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Rayonix4x4(Detector):
    """
    PyFAI Detector instance for the Rayonix 4x4 binning
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000176
        self.pixel_size_um = 176.0
        self.detname = "Rayonix"
        self.calib_shape = (1, 1920, 1920)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 1
        self.asics_shape = (1, 1)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Rayonix5x5(Detector):
    """
    PyFAI Detector instance for the Rayonix 5x5 binning
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000220
        self.pixel_size_um = 220.0
        self.detname = "Rayonix"
        self.calib_shape = (1, 1536, 1536)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 1
        self.asics_shape = (1, 1)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Rayonix6x6(Detector):
    """
    PyFAI Detector instance for the Rayonix 6x6 binning
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000264
        self.pixel_size_um = 264.0
        self.detname = "Rayonix"
        self.calib_shape = (1, 1280, 1280)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 1
        self.asics_shape = (1, 1)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Rayonix8x8(Detector):
    """
    PyFAI Detector instance for the Rayonix 8x8 binning
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000352
        self.pixel_size_um = 352.0
        self.detname = "Rayonix"
        self.calib_shape = (1, 960, 960)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 1
        self.asics_shape = (1, 1)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Rayonix10x10(Detector):
    """
    PyFAI Detector instance for the Rayonix 10x10 binning
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000440
        self.pixel_size_um = 440.0
        self.detname = "Rayonix"
        self.calib_shape = (1, 768, 768)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 1
        self.asics_shape = (1, 1)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class ePix10k2M(Detector):
    """
    PyFAI Detector instance for the ePix10k2M
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.0001
        self.pixel_size_um = 100.0
        self.detname = "ePix10k2M"
        self.calib_shape = (16, 352, 384)
        self.n_modules = self.calib_shape[0]
        self.n_asics = 4
        self.asics_shape = (2, 2)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class ePix10kaQuad(Detector):
    """
    PyFAI Detector instance for the ePix10kaQuad
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.0001
        self.pixel_size_um = 100.0
        self.calib_shape = (4, 352, 384)
        self.detname = "ePix10kaQuad"
        self.n_modules = self.calib_shape[0]
        self.n_asics = 4
        self.asics_shape = (2, 2)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Jungfrau05M(Detector):
    """
    PyFAI Detector instance for the Jungfrau05M
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000075
        self.pixel_size_um = 75.0
        self.calib_shape = (1, 512, 1024)
        self.detname = "Jungfrau05M"
        self.n_modules = self.calib_shape[0]
        self.n_asics = 8
        self.asics_shape = (2, 4)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Jungfrau1M(Detector):
    """
    PyFAI Detector instance for the Jungfrau1M
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000075
        self.pixel_size_um = 75.0
        self.calib_shape = (2, 512, 1024)
        self.detname = "Jungfrau1M"
        self.n_modules = self.calib_shape[0]
        self.n_asics = 8
        self.asics_shape = (2, 4)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Jungfrau4M(Detector):
    """
    PyFAI Detector instance for the Jungfrau4M
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000075
        self.pixel_size_um = 75.0
        self.calib_shape = (8, 512, 1024)
        self.detname = "Jungfrau4M"
        self.n_modules = self.calib_shape[0]
        self.n_asics = 8
        self.asics_shape = (2, 4)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

class Jungfrau16M(Detector):
    """
    PyFAI Detector instance for the Jungfrau16M
    """

    def __init__(
        self,
    ):
        self.pixel_size = 0.000075
        self.pixel_size_um = 75.0
        self.calib_shape = (32, 512, 1024)
        self.detname = "Jungfrau16M"
        self.n_modules = self.calib_shape[0]
        self.n_asics = 8
        self.asics_shape = (2, 4)
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

DETECTOR_REGISTRY = {
    ("Rayonix", (1, 7680, 7680)): Rayonix,
    ("Rayonix", (1, 3840, 3840)): Rayonix2x2,
    ("Rayonix", (1, 2560, 2560)): Rayonix3x3,
    ("Rayonix", (1, 1920, 1920)): Rayonix4x4,
    ("Rayonix", (1, 1536, 1536)): Rayonix5x5,
    ("Rayonix", (1, 1280, 1280)): Rayonix6x6,
    ("Rayonix", (1, 960, 960)):   Rayonix8x8,
    ("Rayonix", (1, 768, 768)):   Rayonix10x10,
    ("ePix10k2M", (16, 352, 384)):   ePix10k2M,
    ("ePix10kaQuad", (4, 352, 384)): ePix10kaQuad,
    ("Jungfrau05M", (1, 512, 1024)):  Jungfrau05M,
    ("Jungfrau1M", (2, 512, 1024)):  Jungfrau1M,
    ("Jungfrau4M", (8, 512, 1024)):  Jungfrau4M,
    ("Jungfrau16M", (32, 512, 1024)): Jungfrau16M,
}

def get_detector(detname: str, shape: tuple) -> Detector:
    """
    Instantiate a PyFAI Detector object based on detector name and shape.

    Parameters
    ----------
    detname : str
        Detector name in psana (e.g., 'jungfrau', 'epix10k2M', etc.)
    shape : tuple
        Unassembled shape
    """
    full_detname = get_pyFAI_detname(detname)
    for (name_pattern, expected_shape), detector_cls in DETECTOR_REGISTRY.items():
        name_matches = full_detname == name_pattern
        shape_matches = shape == expected_shape
        if name_matches and shape_matches:
            Detector.registry[detector_cls.__name__.lower()] = detector_cls
            return detector_cls()

    raise ValueError(
        f"Detector not recognized: detname={detname!r}, shape={shape}. "
        f"Supported detectors: {list(DETECTOR_REGISTRY.keys())}"
    )
from pyFAI.detectors import Detector
from LCLSGeom.utils.detectors import DETECTOR_SPECS, check_detname, get_family

class LCLSDetector(Detector):
    """
    PyFAI Detector instance for an LCLS detector, described by its entry in DETECTOR_SPECS

    Subclasses are generated below, one per entry, so that each keeps its own class name: it is
    the key under which the detector is registered in the PyFAI detector registry.
    """

    def __init__(self):
        self.detname, self.pixel_size_um, self.calib_shape, self.asics_shape = DETECTOR_SPECS[type(self).__name__]
        self.pixel_size = self.pixel_size_um / 1e6
        self.n_modules = self.calib_shape[0]
        self.n_asics = self.asics_shape[0] * self.asics_shape[1]
        self.ss_size = self.calib_shape[1] // self.asics_shape[0]
        self.fs_size = self.calib_shape[2] // self.asics_shape[1]
        super().__init__(pixel1=self.pixel_size, pixel2=self.pixel_size, max_shape=(self.n_modules * self.asics_shape[0] * self.ss_size, self.asics_shape[1] * self.fs_size))

for _name in DETECTOR_SPECS:
    globals()[_name] = type(_name, (LCLSDetector,), {"__doc__": f"PyFAI Detector instance for the {_name}"})

# Keyed by (family, calib shape), the same key check_detname validates against, so the two
# cannot drift: both are derived from DETECTOR_SPECS.
PYFAI_DETECTORS = {
    (family, shape): globals()[name]
    for name, (family, _, shape, _) in DETECTOR_SPECS.items()
}

def get_detector(detname: str, shape: tuple) -> Detector:
    """
    Instantiate the PyFAI Detector matching a detector name and calibration shape

    Parameters
    ----------
    detname : str
        Detector name in psana (e.g., 'jungfrau', 'epix10k2M', etc.)
    shape : tuple
        Unassembled shape as (<number-of-segments>, <rows>, <cols>)

    Returns
    -------
    detector : pyFAI.Detector
        PyFAI Detector instance, also registered in the PyFAI detector registry
    """
    check_detname(detname=detname, shape=shape)
    detector_cls = PYFAI_DETECTORS[(get_family(detname), shape)]
    Detector.registry[detector_cls.__name__.lower()] = detector_cls
    return detector_cls()

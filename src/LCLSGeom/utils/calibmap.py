"""
Calibration group and parameter mappings for LCLS detectors.
"""

import os
from typing import Optional

detname_to_pars = {
    'epix10k2m': 'p0a0,p1a0,p2a0,p3a0,p4a0,p5a0,p6a0,p7a0,p8a0,p9a0,p10a0,p11a0,p12a0,p13a0,p14a0,p15a0',
    'epix10kaquad': 'p0a0,p1a0,p2a0,p3a0',
    'jungfrau05m': 'p0a0',
    'jungfrau1m': 'p0a0,p1a0',
    'jungfrau4m': 'p0a0,p1a0,p2a0,p3a0,p4a0,p5a0,p6a0,p7a0',
    'jungfrau16m': 'p0a0,p1a0,p2a0,p3a0,p4a0,p5a0,p6a0,p7a0,p8a0,p9a0,p10a0,p11a0,p12a0,p13a0,p14a0,p15a0',
}

detname_to_groups = {
    'epix10k2m': {
                  'group_p0': 'p0a0,p0a1,p0a2,p0a3', 'group_p1': 'p1a0,p1a1,p1a2,p1a3',
                  'group_p2': 'p2a0,p2a1,p2a2,p2a3', 'group_p3': 'p3a0,p3a1,p3a2,p3a3',
                  'group_p4': 'p4a0,p4a1,p4a2,p4a3', 'group_p5': 'p5a0,p5a1,p5a2,p5a3',
                  'group_p6': 'p6a0,p6a1,p6a2,p6a3', 'group_p7': 'p7a0,p7a1,p7a2,p7a3',
                  'group_p8': 'p8a0,p8a1,p8a2,p8a3', 'group_p9': 'p9a0,p9a1,p9a2,p9a3',
                  'group_p10': 'p10a0,p10a1,p10a2,p10a3', 'group_p11': 'p11a0,p11a1,p11a2,p11a3',
                  'group_p12': 'p12a0,p12a1,p12a2,p12a3', 'group_p13': 'p13a0,p13a1,p13a2,p13a3',
                  'group_p14': 'p14a0,p14a1,p14a2,p14a3', 'group_p15': 'p15a0,p15a1,p15a2,p15a3',
                  },
    'epix10kaquad': {
                     'group_p0': 'p0a0,p0a1,p0a2,p0a3', 'group_p1': 'p1a0,p1a1,p1a2,p1a3',
                     'group_p2': 'p2a0,p2a1,p2a2,p2a3', 'group_p3': 'p3a0,p3a1,p3a2,p3a3',
                    },
    'jungfrau05m': {
                    'group_p0': 'p0a0,p0a1,p0a2,p0a3,p0a4,p0a5,p0a6,p0a7'
                   },
    'jungfrau1m': {
                   'group_p0': 'p0a0,p0a1,p0a2,p0a3,p0a4,p0a5,p0a6,p0a7', 
                   'group_p1': 'p1a0,p1a1,p1a2,p1a3,p1a4,p1a5,p1a6,p1a7',
                  },
    'jungfrau4m': {
                   'group_p0': 'p0a0,p0a1,p0a2,p0a3,p0a4,p0a5,p0a6,p0a7',
                   'group_p1': 'p1a0,p1a1,p1a2,p1a3,p1a4,p1a5,p1a6,p1a7',
                   'group_p2': 'p2a0,p2a1,p2a2,p2a3,p2a4,p2a5,p2a6,p2a7',
                   'group_p3': 'p3a0,p3a1,p3a2,p3a3,p3a4,p3a5,p3a6,p3a7',
                   'group_p4': 'p4a0,p4a1,p4a2,p4a3,p4a4,p4a5,p4a6,p4a7',
                   'group_p5': 'p5a0,p5a1,p5a2,p5a3,p5a4,p5a5,p5a6,p5a7',
                   'group_p6': 'p6a0,p6a1,p6a2,p6a3,p6a4,p6a5,p6a6,p6a7',
                   'group_p7': 'p7a0,p7a1,p7a2,p7a3,p7a4,p7a5,p7a6,p7a7',
                  },
    'jungfrau16m': {
                    'group_p0': 'p0a0,p0a1,p0a2,p0a3,p0a4,p0a5,p0a6,p0a7',
                    'group_p1': 'p1a0,p1a1,p1a2,p1a3,p1a4,p1a5,p1a6,p1a7',
                    'group_p2': 'p2a0,p2a1,p2a2,p2a3,p2a4,p2a5,p2a6,p2a7',
                    'group_p3': 'p3a0,p3a1,p3a2,p3a3,p3a4,p3a5,p3a6,p3a7',
                    'group_p4': 'p4a0,p4a1,p4a2,p4a3,p4a4,p4a5,p4a6,p4a7',
                    'group_p5': 'p5a0,p5a1,p5a2,p5a3,p5a4,p5a5,p5a6,p5a7',
                    'group_p6': 'p6a0,p6a1,p6a2,p6a3,p6a4,p6a5,p6a6,p6a7',
                    'group_p7': 'p7a0,p7a1,p7a2,p7a3,p7a4,p7a5,p7a6,p7a7',
                    'group_p8': 'p8a0,p8a1,p8a2,p8a3,p8a4,p8a5,p8a6,p8a7',
                    'group_p9': 'p9a0,p9a1,p9a2,p9a3,p9a4,p9a5,p9a6,p9a7',
                    'group_p10': 'p10a0,p10a1,p10a2,p10a3,p10a4,p10a5,p10a6,p10a7',
                    'group_p11': 'p11a0,p11a1,p11a2,p11a3,p11a4,p11a5,p11a6,p11a7',
                    'group_p12': 'p12a0,p12a1,p12a2,p12a3,p12a4,p12a5,p12a6,p12a7',
                    'group_p13': 'p13a0,p13a1,p13a2,p13a3,p13a4,p13a5,p13a6,p13a7',
                    'group_p14': 'p14a0,p14a1,p14a2,p14a3,p14a4,p14a5,p14a6,p14a7',
                    'group_p15': 'p15a0,p15a1,p15a2,p15a3,p15a4,p15a5,p15a6,p15a7',
                    'group_p16': 'p16a0,p16a1,p16a2,p16a3,p16a4,p16a5,p16a6,p16a7',
                    'group_p17': 'p17a0,p17a1,p17a2,p17a3,p17a4,p17a5,p17a6,p17a7',
                    'group_p18': 'p18a0,p18a1,p18a2,p18a3,p18a4,p18a5,p18a6,p18a7',
                    'group_p19': 'p19a0,p19a1,p19a2,p19a3,p19a4,p19a5,p19a6,p19a7',
                    'group_p20': 'p20a0,p20a1,p20a2,p20a3,p20a4,p20a5,p20a6,p20a7',
                    'group_p21': 'p21a0,p21a1,p21a2,p21a3,p21a4,p21a5,p21a6,p21a7',
                    'group_p22': 'p22a0,p22a1,p22a2,p22a3,p22a4,p22a5,p22a6,p22a7',
                    'group_p23': 'p23a0,p23a1,p23a2,p23a3,p23a4,p23a5,p23a6,p23a7',
                    'group_p24': 'p24a0,p24a1,p24a2,p24a3,p24a4,p24a5,p24a6,p24a7',
                    'group_p25': 'p25a0,p25a1,p25a2,p25a3,p25a4,p25a5,p25a6,p25a7',
                    'group_p26': 'p26a0,p26a1,p26a2,p26a3,p26a4,p26a5,p26a6,p26a7',
                    'group_p27': 'p27a0,p27a1,p27a2,p27a3,p27a4,p27a5,p27a6,p27a7',
                    'group_p28': 'p28a0,p28a1,p28a2,p28a3,p28a4,p28a5,p28a6,p28a7',
                    'group_p29': 'p29a0,p29a1,p29a2,p29a3,p29a4,p29a5,p29a6,p29a7',
                    'group_p30': 'p30a0,p30a1,p30a2,p30a3,p30a4,p30a5,p30a6,p30a7',
                    'group_p31': 'p31a0,p31a1,p31a2,p31a3,p31a4,p31a5,p31a6,p31a7',
                   },
}

detname_to_quadrants = {
    'epix10k2m': {
                  'group_q0': 'p0,p1,p2,p3', 
                  'group_q1': 'p4,p5,p6,p7',
                  'group_q2': 'p8,p9,p10,p11', 
                  'group_q3': 'p12,p13,p14,p15',
                  'group_all': 'q0,q1,q2,q3',
                 },
    'epix10kaquad': None,
    'jungfrau05m': None,
    'jungfrau1m': None,
    'jungfrau4m': {
                   'group_q0': 'p0,p1',
                   'group_q1': 'p2,p3',
                   'group_q2': 'p4,p5',
                   'group_q3': 'p6,p7',
                   'group_all': 'q0,q1,q2,q3',
                  },
    'jungfrau16m': {
                    'group_q0': 'p0,p1,p2,p3,p4,p5,p6,p7',
                    'group_q1': 'p8,p9,p10,p11,p12,p13,p14,p15',
                    'group_q2': 'p16,p17,p18,p19,p20,p21,p22,p23',
                    'group_q3': 'p24,p25,p26,p27,p28,p29,p30,p31',
                    'group_all': 'q0,q1,q2,q3',
                   },
}

calib_groups = (
    "UNDEFINED",
    "Camera::CalibV1",
    "Jungfrau::CalibV1",
    "Jungfrau::CalibV1",
    "Jungfrau::CalibV1",
    "Jungfrau::CalibV1",
    "Epix10ka2M::CalibV1",
    "Epix10kaQuad::CalibV1",
    "Epix10kaQuad::CalibV1",
    "Epix10kaQuad::CalibV1",
    "Epix10kaQuad::CalibV1",
)

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

full_detnames = (
    "UNDEFINED",
    "rayonix",
    "jungfrau05M",
    "jungfrau1M",
    "jungfrau4M",
    "jungfrau16M",
    "epix10k2M",
    "epix10kaQuad0",
    "epix10kaQuad1",
    "epix10kaQuad2",
    "epix10kaQuad3",
)

calib_detnames = (
    "UNDEFINED",
    "Rayonix",
    "Jungfrau",
    "Jungfrau",
    "Jungfrau",
    "Jungfrau",
    "Epix10ka2M",
    "Epix10kaQuad.0",
    "Epix10kaQuad.1",
    "Epix10kaQuad.2",
    "Epix10kaQuad.3",
)

pyFAI_detnames = (
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

hutches = (
    "UNDEFINED",
    "XPP",
    "XCS",
    "CXI",
    "MEC",
    "MFX",
)

stations = (
    "UNDEFINED",
    "XppEndstation.0",
    "XcsEndstation.0",
    "CxiDs1.0",
    "MecTargetChamber.0",
    "MfxEndstation.0",
)

psana_detnames_lower = tuple(name.lower() for name in psana_detnames)

psana_to_calib_detname = dict(zip(psana_detnames_lower, calib_detnames))

psana_to_full_detname = dict(zip(psana_detnames_lower, full_detnames))

psana_to_pyFAI_detname = dict(zip(psana_detnames_lower, pyFAI_detnames))

det_to_group = dict(zip(calib_detnames, calib_groups))

hutch_to_station = dict(zip(hutches, stations))


def group_from_detname(detname: str) -> str:
    """Retrieve the group string from the detector type."""
    detname = psana_to_calib_detname.get(detname.lower(), "UNDEFINED")
    if detname == "UNDEFINED":
        raise ValueError(f"Unknown detector type: {detname}")
    group = det_to_group.get(detname)
    return group


def source_from_detname(detname: str, hutch: str) -> str:
    """Retrieve the source string from the detector name and hutch."""
    station = hutch_to_station.get(hutch.upper(), "UNDEFINED")
    if station == "UNDEFINED":
        raise ValueError(f"Unknown hutch: {hutch}")
    detname = psana_to_calib_detname.get(detname.lower(), "UNDEFINED")
    if detname == "UNDEFINED":
        raise ValueError(f"Unknown detector type: {detname}")
    if "Epix10kaQuad" in detname:
        return f"{station}:{detname}"
    return f"{station}:{detname}.0"


def select_calib_file(calib_dir: str, run: int) -> Optional[str]:
    """Select the calibration file from the calibration directory and run number."""
    if not os.path.exists(calib_dir):
        return

    fnames = os.listdir(calib_dir)
    files = [os.path.join(calib_dir, fname) for fname in fnames]

    run_max = 9999
    run_files = []
    for file in files:
        f = os.path.basename(file)
        if f == "HISTORY":
            continue
        if os.path.splitext(f)[1] != ".data":
            continue
        basename = os.path.splitext(f)[0]
        fields = basename.split("-")
        begin, end = fields

        if begin.isdigit():
            begin_int = int(begin)
            if begin_int >= run_max:
                raise ValueError(
                    f"Begin run number {run} is too high for calibration directory {calib_dir}"
                )

        if end.isdigit():
            end_int = int(end)
            if end_int >= run_max:
                raise ValueError(
                    f"End run number {run} is too high for calibration directory {calib_dir}"
                )
        elif end == "end":
            end_int = run_max

        run_files.append((begin_int, end_int, file))
    run_files.sort(key=lambda x: int(x[0]))

    for run_file in run_files[::-1]:
        if run_file[0] <= run <= run_file[1]:
            return run_file[2]

    return os.path.join(calib_dir, "0-end.data")


def get_full_detname(detname: str) -> str:
    """Start from psana detector name and return the full detector name."""
    full_detname = psana_to_full_detname.get(detname.lower(), "UNDEFINED")
    if full_detname == "UNDEFINED":
        raise ValueError(f"Unknown detector type: {detname}")
    return full_detname

def get_pyFAI_detname(detname: str) -> str:
    """Start from psana detector name and return the pyFAI detector name."""
    pyFAI_detname = psana_to_pyFAI_detname.get(detname.lower(), "UNDEFINED")
    if pyFAI_detname == "UNDEFINED":
        raise ValueError(f"Unknown detector type: {detname}")
    return pyFAI_detname
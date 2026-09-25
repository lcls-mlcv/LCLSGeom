# LCLSGeom
LCLSwitching Geometry tool - helper classes and functions for managing and converting geometry files for LCLS-I and LCLS-II detectors.
This repository aims at providing an easy-to-use geometry management and conversion tool for LCLS scripts. The source code mainly relies on [PSCalib](https://github.com/lcls-psana/PSCalib/tree/master) for LCLS-I detectors and [psana](https://github.com/slac-lcls/lcls2/tree/master) for LCLS-II detectors, both repositeries designed by Mikhail Dubrovin.

# Table of Contents
- [Understanding Geometry Frameworks](#understanding-geometry-frameworks)
  - [Psana and CrystFEL](#psana-and-crystfel)
  - [PyFAI](#pyfai)
  - [LCLSGeom Motivation](#lclsgeom-motivation)
- [User Guide](#user-guide)
  - [Installation](#installation)
  - [Package Modules](#package-modules)
  - [Command-line Interface](#command-line-interface)
- [Organisation](#organization)

## Understanding Geometry Frameworks 

At LCLS, many geometry formats coexist. Some workflows rely on `CrystFEL` to process data, some others rely on `PyFAI`. This is without mentioning LCLS internal geometry framework which will be referred as `psana` in this package. This project aims to centralize the conversion scripts between those formats and to provide helper functions to interact with the geometry calibration LCLS instances.

### Psana & CrystFEL

Both psana and CrystFEL geometry formats describe detectors using a hierarchical 3D structure. In this representation, the detector geometry is defined through multiple coordinate frames arranged in a hierarchy:
```
Interaction Point (IP)
        ↓
  Detector Center
        ↓
      Panels
        ↓
      Asics
```

Each level defines the position and orientation of its children relative to its own coordinate frame. Leaf objects (panels or ASICs) are therefore described by their translation and rotation with respect to the detector center, which itself is positioned relative to the interaction point.

This approach is particularly well suited for modern pixel detectors used at LCLS, which are typically segmented detectors composed of multiple panels arranged in 3D space. The hierarchical structure allows accurate modeling of panel offsets, tilts, and metrology corrections.

For example, here is the structure of a `psana` .data file:

```
# TITLE      Geometry parameters of ePix10k2M
# DATE_TIME  2025-01-01T00:00:00 PST
# METROLOGY  no metrology available
# AUTHOR     dubrovin
# EXPERIMENT any
# DETECTOR   epix10k2M
# CALIB_TYPE geometry
# COMMENT:01 Manually created as an example for the 16-segment ePix10k2M detector

# PARAM:01 PARENT     - name and version of the parent object
# PARAM:02 PARENT_IND - index of the parent object
# PARAM:03 OBJECT     - name and version of the object
# PARAM:04 OBJECT_IND - index of the new object
# PARAM:05 X0         - x-coordinate [um] of the object origin in the parent frame
# PARAM:06 Y0         - y-coordinate [um] of the object origin in the parent frame
# PARAM:07 Z0         - z-coordinate [um] of the object origin in the parent frame
# PARAM:08 ROT_Z      - object design rotation angle [deg] around Z axis of the parent frame
# PARAM:09 ROT_Y      - object design rotation angle [deg] around Y axis of the parent frame
# PARAM:10 ROT_X      - object design rotation angle [deg] around X axis of the parent frame
# PARAM:11 TILT_Z     - object tilt angle [deg] around Z axis of the parent frame
# PARAM:12 TILT_Y     - object tilt angle [deg] around Y axis of the parent frame
# PARAM:13 TILT_X     - object tilt angle [deg] around X axis of the parent frame

# HDR PARENT IND     OBJECT IND    X0[um]   Y0[um]   Z0[um]   ROT-Z  ROT-Y  ROT-X     TILT-Z    TILT-Y    TILT-X

      CAMERA  0  EPIX10KA:V1  0    -58583    24251        5     270      0      0   -1.19826   0.11768  -0.03984 # This line sets where is Panel 0 relative to CAMERA
      CAMERA  0  EPIX10KA:V1  1    -57864    64862       66     270      0      0   -1.16980  -0.22605   0.07420 # ...
      CAMERA  0  EPIX10KA:V1  2    -16705    23716      -87     270      0      0   -0.88179  -0.04979  -0.06790 # ...
      CAMERA  0  EPIX10KA:V1  3    -16091    64293        1     270      0      0   -1.12086  -0.09023  -0.10545 # ...
      CAMERA  0  EPIX10KA:V1  4     24315    58774       -4     180      0      0   -0.94978  -0.01083  -0.15229 # ...
      CAMERA  0  EPIX10KA:V1  5     64818    58098       -5     180      0      0   -0.94860  -0.02310  -0.14836 # ...
      CAMERA  0  EPIX10KA:V1  6     23737    16960       35     180      0      0   -0.81313   0.10535  -0.22019 # ...
      CAMERA  0  EPIX10KA:V1  7     64176    16289       15     180      0      0   -0.89935  -0.13787   0.18512 # ...
      CAMERA  0  EPIX10KA:V1  8     58761   -24236      113      90      0      0   -1.18462  -0.05487  -0.00937 # This line sets where is Panel i relative to CAMERA
      CAMERA  0  EPIX10KA:V1  9     58014   -65006      -75      90      0      0   -1.09273   0.06857  -0.19293 # ...
      CAMERA  0  EPIX10KA:V1 10     16766   -23772     -150      90      0      0   -0.91635  -0.15069   0.04136 # ...
      CAMERA  0  EPIX10KA:V1 11     16175   -64142       59      90      0      0   -0.92756  -0.14366  -0.04918 # ...
      CAMERA  0  EPIX10KA:V1 12    -24466   -58939       17       0      0      0   -1.24074  -0.00144  -0.03904 # ...
      CAMERA  0  EPIX10KA:V1 13    -65201   -58201       32       0      0      0   -0.90080   0.02960   0.23655 # ...
      CAMERA  0  EPIX10KA:V1 14    -23760   -16774      -49       0      0      0   -0.66263  -0.27791  -0.04138 # ...
      CAMERA  0  EPIX10KA:V1 15    -64091   -16171       24       0      0      0   -1.09707   0.02238  -0.01484 # ...
          IP  0       CAMERA  0         0    -1000   100000      90      0      0    0.00000   0.00000   0.00000 # This line sets where CAMERA is relative to IP
```

As you can notice, this example .data file specifies an offset in the fast dimension of -0.001 m, a detector-sample distance of 0.1 m and that the detector is rotated by 90 degrees around the detector's normal.

### PyFAI

Most PyFAI workflows, however, assume a single monolithic 2D detector with a regular pixel grid. In this model, pixels are indexed directly in a 2D grid and their positions are inferred from the detector geometry parameters. To be precise, 6 geometry parameters  (sample-PONI distance `dist`, PONI (for point of normal incidence) coordinates `poni1` and `poni2`, rotations around 1st, 2nd and 3rd axis `rot1`, `rot2` and `rot3`) *define* the geometry.

For segmented detectors, this means that raw detector data must first be assembled into a single 2D image before PyFAI can operate on it. Assembly is the process of projecting the individual detector panels from their 3D positions into a continuous 2D pixel layout.

While this approach works, it introduces an additional preprocessing step that must be repeated each time an image is processed, which can become computationally inefficient when dealing with large datasets. Furthermore, the assembly step implicitly projects the detector geometry onto a single plane. For detectors that are not perfectly planar, this projection forces a non-planar geometry into a planar representation, which can introduce small geometric distortions and artefacts in the reconstructed image and potentially propagate into downstream analysis.

### LCLSGeom Motivation

LCLSGeom aims to bridge the gap between these two geometry paradigms.

Instead of assembling the detector image at runtime, LCLSGeom constructs a virtual 2D detector representation compatible with PyFAI. This is achieved by creating a fake 2D detector where the panel dimensions are expanded along the slow-scan direction so that all panel pixels can be mapped into a single continuous 2D indexing scheme.

The true physical pixel positions are then encoded through the detector metrology information, which is passed to PyFAI. This allows PyFAI to correctly interpret the spatial layout of the pixels while still operating on a single 2D detector structure.

The metrology information can originate from any .data geometry file. As a result, *the PONI parameters are defined relatively to the metrology used as a reference*. If the input .data file presents axis offsets (such as a detector position which is not centered on the detector center for example), the PyFAI PONI parameters will then be defined relative to those offsets. Hence, implicetely, LCLSGeom will zero-out any offsets present in the metrology, unless specified the contrary. This is done by setting the parameter `image_frame=False` to any conversion scripts involving `PyFAI`.

By encoding the detector layout in this way, LCLSGeom allows PyFAI to operate directly on detector data without requiring explicit image assembly, improving efficiency while preserving accurate detector geometry.

## User Guide

### Installation

LCLSGeom is available in both the LCLS-I (`psana1`) and LCLS-II (`psana2`) environments and does not require a separate installation in those contexts. To activate the appropriate environment on S3DF:
```bash
# LCLS-I (psana1)
source /sdf/group/lcls/ds/ana/sw/conda1/manage/bin/psconda.sh

# LCLS-II (psana2)
source /sdf/group/lcls/ds/ana/sw/conda2/manage/bin/psconda.sh
```

For use outside of these environments, the package can be installed directly from the repository:
```bash
git clone https://github.com/slac-lcls/LCLSGeom.git
cd LCLSGeom
pip install -e .
```

> **Nota Bene:** LCLSGeom depends on `psana` and `PSCalib`, which are not available on PyPI and cannot be installed via `pip`.
> Reading a `psana` geometry file requires `GeometryAccess`, so **all** conversions — including `lcls-convert` — need a `psana1` or `psana2` environment.
> A standalone `pip` installation is therefore of limited use on its own.
> Conversions to `dxtbx` additionally require `dxtbx`, which is conda-only; it is imported lazily, so the rest of the package works without it. `pyFAI` is likewise imported lazily, so the `dxtbx` and `CrystFEL` converters do not need it.
> The `dxtbx` converters are available in the LCLS `cctbx` builds, which provide `dxtbx` and `psana` together:
> ```bash
> source /sdf/group/lcls/ds/tools/cctbx/psana2_setup.sh   # or psana1_setup.sh
> ```

### Package Modules

LCLSGeom consists of two modules:
- _manager_: Handles relations with calibration databases (for LCLS-II detectors) or directories (for LCLS-I detectors).
- _converter_: Provides quick conversions from one format to another. Available formats are `psana`, `CrystFEL`, `pyFAI` and `dxtbx`.

#### From Psana Geometry to PyFAI detectors

The _PsanaToPyFAI_ class is able to parse a `psana` geometry file and converts the geometry information in a `pyFAI` detector object that can be used to perform any downstream `pyFAI` task.

##### Example Usage

**Import required modules**
```python
from LCLSGeom.manager import get_geometry
from LCLSGeom.converter import PsanaToPyFAI
```

LCLS-I experiment
```python
exp = mfxx49820
run = 8
detname = "epix10k2M"

geometry_file = get_geometry(detname, exp, run) # <-- Will fetch geometry file stored in calibration directory or on LCLSGeom/templates if not found
epix10k2M_detector = PsanaToPyFAI.convert(geometry_file, detname) # <-- Implicitly `image_frame=True` zero-out offsets to set up coordinates in image frame
x, y, z = epix10k2M_detector.calc_cartesian_positions()
x.mean()
```
Output
```python
0.0
```

LCLS-II experiment
```python
exp = "mfx101591226"
run = 20
detname = "jungfrau"

geometry_file = get_geometry(detname, exp, run) # <-- Will fetch geometry file stored in slac-lcls/lcls2 or on LCLSGeom/templates if not found
jungfrau_detector = PsanaToPyFAI.convert(geometry_file, detname) # <-- Implicitly `image_frame=True` zero-out offsets to set up coordinates in image frame
x, y, z = jungfrau_detector.calc_cartesian_positions()
x.mean()
```
Output
```python
0.0
```

> **Nota Bene:** By default, PsanaToPyFAI.convert removes any axis offsets present in the metrology .data file. 
> This ensures that the resulting PyFAI geometry is centered at (x, y) = (0, 0) and z = 0. 
> In other words, the detector origin is placed at the image center, rather than at the interaction point (IP).
> If needed, this behavior can be disabled by passing `image_frame=False`. In that case, the PyFAI geometry is defined in the laboratory coordinate frame.

###### Bonus Example Usage

If we convert the original .data example with an offset in the slow dimension, here is what we can get by passing `image_frame=False`.
 ```python
exp = mfxx49820
run = 8
detname = "epix10k2M"

geometry_file = get_geometry(detname) # <-- Will fetch default geometry file stored in LCLSGeom/templates
epix10k2M_detector = PsanaToPyFAI.convert(geometry_file, detname, image_frame=False) # <-- Do not zero-out offsets to set up coordinates in laboratory frame
x, y, z = jungfrau_detector.calc_cartesian_positions()
x.mean()
```
Output
```python
-0.001
```

#### From PyFAI detectors back to Psana

The _PyFAIToPsana_ class is able to write a `psana` geometry file using a `PyFAI` detector object and a corresponding .poni file containing the PONI parameters necessary to define the entire detector geometry.
This class is designed to update geometries after a pyFAI geometry optimization routine was called.

##### Example Usage

**Import required modules**
```python
from LCLSGeom.converter import PsanaToPyFAI, PyFAIToPsana
```

LCLS-I or LCLS-II experiment
```python
exp = "mfxx49820"
run = 8
detname = "epix10k2M"

geometry_file = get_geometry(detname) # <-- Will fetch default geometry file
epix10k2M_detector = PsanaToPyFAI.convert(geometry_file, detname) # <-- Implicitly `image_frame=True` zero-out offsets to set up coordinates in image frame
# ----- Analysis ------
# ...
# ------- Done! -------
poni_file = "path/to/file.poni" # <-- This is the result of the geometry calibration using PyFAI
out_file = "path/to/desired/folder/<run>-end.data" # <-- New geometry
PyFAIToPsana.convert(poni_file, epix10k2M_detector, out_file) # <-- Implicitly `image_frame=True` will write .data file in image frame
```

By default, the new geometry will be defined in the image frame, i.e, no 90° degree rotation around z-axis will be applied to match `psana` conventions.
This can be disabled by passing `image_frame=False` to output the geometry file in the `psana` coordinate system.

```python
# ----- Analysis ------
# ...
# ------- Done! -------
poni_file = "path/to/file.poni" # <-- This is the result of the geometry calibration using PyFAI
out_file = "path/to/desired/folder/<run>-end.data" # <-- New geometry
PyFAIToPsana.convert(poni_file, epix10k2M_detector, out_file, image_frame=False) # <-- This will add the 90 degree rotation at the IP-detector level
```

#### Switching between Psana and CrystFEL

The _PsanaToCrystFEL_ class is able to parse a `psana` geometry file and converts the geometry information into a `CrystFEL` geometry file.
The _CrystFELToPsana_ class is able to parse a `CrystFEL` geometry file and converts it back to a `psana` geometry file.

##### Example Usage
```python
from LCLSGeom.manager import get_geometry
from LCLSGeom.converter import PsanaToCrystFEL

# LCLS-I or LCLS-II experiment
exp = mfxx49820
run = 8
detname = "epix10k2M"

data_file = get_geometry(detname, exp, run)
geom_file = "path/to/desired/file.geom"
PsanaToCrystFEL.convert(data_file, detname, geom_file)
data_file_2 = "path/to/desired/file.data"
CrystFELToPsana.convert(geom_file, detname, data_file_2)
```

#### From Psana Geometry to dxtbx detectors

Two converters produce `dxtbx` `Detector` models, for the two incompatible conventions used
downstream. Both work in millimetres, since that is what `dxtbx` expects, whereas `psana`
reports microns.

| | _PsanaToDXTBX_ | _PsanaToCCTBX_ |
|---|---|---|
| Intended consumer | `simtbx.nanoBragg` | `DIALS` / `cctbx.xfel` |
| Panels | one per ASIC | as `cctbx` does: padded module-level for Jungfrau, per-ASIC otherwise |
| Jungfrau4M panel count / size | 64 × (256, 256) | 8 × (1030, 514) |
| Structure | flat, absolute origins | hierarchical, origins local to each module group |
| Detector placement | zeroed | preserved from the input file |
| Raw data | usable directly | needs `serialtbx.detector.jungfrau.correct_panel` |
| `convert(...)` | `(in_file, detname)` → `Detector` | `(in_file, detname)` → `Detector` |
| `convert_to_file(...)` | `(in_file, detname, out_file)` → writes `.expt` | `(in_file, detname, out_file)` → writes `.expt` |

_PsanaToDXTBX_ emits one panel per ASIC, matching the `psana` calibration array one-to-one, so
raw frames need no padding. All detector placement is zeroed — the camera offsets, the 90°
rotation of the `psana` convention, and the sample-detector distance — so panel origins hold
pure metrology and placement is supplied downstream as geometric parameters. Because the total
distance is then zero, `get_beam_centre()` and `get_ray_intersection()` cannot be evaluated on
a planar detector.

_PsanaToCCTBX_ reproduces what `cctbx.xfel` writes, so its output is directly comparable with
refined `.expt` files. Nothing is zeroed. Jungfrau segments are emitted as single panels padded
for the double-width pixels at the inter-ASIC boundaries (512×1024 → 514×1030), which means raw
frames must be passed through `serialtbx.detector.jungfrau.correct_panel` first.

##### Example Usage
```python
from LCLSGeom.manager import get_geometry
from LCLSGeom.converter import PsanaToDXTBX, PsanaToCCTBX

exp = "mfx101591226"
run = 20
detname = "jungfrau"

data_file = get_geometry(detname, exp, run)

# Flat, per-ASIC panels for simulation: returns a dxtbx Detector
detector = PsanaToDXTBX.convert(data_file, detname)

# ...or write the same geometry to an .expt file
PsanaToDXTBX.convert_to_file(data_file, detname, "path/to/desired/file.expt")

# cctbx-convention geometry for DIALS: returns a dxtbx Detector
detector = PsanaToCCTBX.convert(data_file, detname)

# ...or write the same geometry to an .expt file
PsanaToCCTBX.convert_to_file(data_file, detname, "path/to/desired/file.expt")
```

### Command-line Interface

LCLSGeom provides two command-line tools.

**`lcls-convert`** converts geometry files between `psana`, `CrystFEL` and `dxtbx` formats. The conversion direction is inferred automatically from the file extensions (`.data` ↔ `.geom`, `.data` → `.json`/`.expt`).
```bash
# psana → CrystFEL
lcls-convert -i input.data -d epix10k2M -o output.geom

# CrystFEL → psana
lcls-convert -i input.geom -d jungfrau -o output.data

# psana → dxtbx, flat per-ASIC panels (simtbx.nanoBragg)
lcls-convert -i input.data -d jungfrau4M -o output.expt

# psana → dxtbx, cctbx-convention hierarchy (DIALS / cctbx.xfel)
lcls-convert -i input.data -d jungfrau4M -o output.expt --cctbx
```

| Argument | Description |
|---|---|
| `--cctbx` | For `dxtbx` output, emit the `cctbx` module-level hierarchical geometry instead of flat per-ASIC panels |

**`lcls-push-db`** pushes a geometry file to the LCLS calibration database. This requires a `psana2` environment.
```bash
lcls-push-db -e mfxx49820 -r 8 -d epix10k2M -g path/to/geometry.data
```

| Argument | Description |
|---|---|
| `-e` | Experiment name |
| `-r` | Run number |
| `-d` | Detector name |
| `-g` | Path to the geometry file |

## Organization

The package structure is organized as follows:

```
LCLSGeom/
├── src/LCLSGeom/
│ ├── __init__.py
│ ├── converter.py            # The converter classes: the public API
│ ├── manager.py              # Geometry file lookup: calib DB, calib directory, templates
│ ├── detector.py             # PyFAI Detector subclasses (the only module importing pyFAI eagerly)
│ ├── utils/
│ │ ├── __init__.py
│ │ ├── detectors.py          # Physical detector facts, families, and legal calibration shapes
│ │ ├── pscalib_utils.py      # psana names, the pscalib geometry tree, and the .data format
│ │ ├── crystfel_utils.py     # The CrystFEL .geom format and its millepede groups
│ │ ├── pyfai_utils.py        # PONI loading and correction, pyFAI frame transforms
│ │ ├── dxtbx_utils.py        # dxtbx/cctbx detector model helpers
│ │ └── geomcalc.py           # Geometry helpers shared by every converter
│ ├── scripts/
│ │ ├── __init__.py
│ │ ├── convert.py            # Command-line Script for Conversion
│ │ └── push_to_database.py   # Command-line Script for Pushing to Database
│ ├── templates/    # Template files for Defaulting Geometries
│ │ ├── geometry-def-epix10k2M.data
│ │ ├── geometry-def-epix10kaQuad0.data
│ │ ├── geometry-def-epix10kaQuad1.data
│ │ ├── geometry-def-epix10kaQuad2.data
│ │ ├── geometry-def-epix10kaQuad3.data
│ │ ├── geometry-def-jungfrau1M.data
│ │ ├── geometry-def-jungfrau4M.data
│ │ ├── geometry-def-jungfrau05M.data
│ │ ├── geometry-def-jungfrau16M.data
│ │ └── geometry-def-rayonix.data
├── tests/
│ ├── test_pyfai.py           # Needs a psana environment with pyFAI
│ └── test_dxtbx.py           # Needs a cctbx environment with psana and dxtbx
├── LICENSE
├── README.md
└── pyproject.toml # Dependencies and requirements
```

Everything under `utils/` is split by the software it speaks to: one module per output format,
each owning its own detector-name convention, plus `detectors.py` for the facts that belong to
the hardware rather than to any one program. No module under `utils/` imports `pyFAI`, `dxtbx`
or `psana` at module level, so the package is importable in a `psana` environment without
`dxtbx` and in a `cctbx` environment without `pyFAI`.
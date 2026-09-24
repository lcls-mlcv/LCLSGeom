"""
Converter classes to convert between psana .data geometry files, CrystFEL .geom geometry files,
PyFAI Detector models, and dxtbx Detector models
"""

try:
    from psana.pscalib.geometry.GeometryAccess import GeometryAccess  # psana2
except ImportError:
    from PSCalib.GeometryAccess import GeometryAccess  # psana1

import os
import numpy as np
from math import atan2, degrees
import pyFAI
from LCLSGeom.manager import get_default_geometry
from LCLSGeom.detector import get_detector
from LCLSGeom.utils.calibmap import detname_to_pars, psana_to_pyFAI_detname
from LCLSGeom.utils.formatting import str_is_segment_and_asic, sfields_to_xyz_vector, header_psana, header_crystfel, tailer_crystfel
from LCLSGeom.utils.frame import image_to_pyfai, pyfai_to_image
from LCLSGeom.utils.geomcalc import correct_geometry, angle_and_tilt, tilt_xy
pyFAI.use_opencl = False

def check_detname(detname, shape):
    """
    Check that a detector name is consistent with the shape of its geometry file

    The dxtbx converters do not need a pyFAI Detector, but the registry is still the
    authoritative map from a psana detector name to its expected calibration shape, so it is
    used here to reject a name that disagrees with the file it was passed.

    Parameters
    ----------
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
    shape : tuple
        Shape of the geometry file as (<number-of-segments>, <rows>, <cols>)
    """
    get_detector(detname=detname, shape=shape)

def write_experiment_list(detector, out_file):
    """
    Write a dxtbx Detector to a geometry file as a one-experiment ExperimentList

    dxtbx and DIALS read .expt files as an ExperimentList; a bare detector dictionary is not
    loadable, so the detector is wrapped in a single Experiment before being written.

    Parameters
    ----------
    detector : dxtbx.model.Detector
        Detector to serialize
    out_file : str
        Path to the output dxtbx .expt/.json geometry file
    """
    from dxtbx.model import Experiment, ExperimentList

    experiments = ExperimentList([Experiment(detector=detector)])
    experiments.as_file(out_file)
    print(f"Wrote dxtbx geometry file {out_file}.", flush=True)

class PsanaToCrystFEL:
    """
    Class to convert psana .data geometry files to CrystFEL .geom geometry files

    Parameters
    ----------
    in_file: str
        Path to the input psana .data geometry file
    """

    def __init__(self, in_file):
        self.geo = GeometryAccess(path=in_file, pbits=0, use_wide_pix_center=False)

    @classmethod
    def convert(cls, in_file, detname, out_file):
        """
        Class method to convert psana .data geometry file to CrystFEL .geom geometry file

        Parameters
        ----------
        in_file : str
            Path to the input psana .data geometry file
        detname : str
            Detector name in psana (e.g., 'jungfrau', 'epix10k2M', etc.)
        out_file : str
            Path to the output CrystFEL .geom file
        """
        converter = cls(in_file)
        converter.convert_data_to_geom(detname=detname, out_file=out_file)

    def convert_data_to_geom(self, detname, out_file):
        """
        Write a CrystFEL .geom file from a psana .data file

        Parameters
        ----------
        detname : str
            Detector name in psana (e.g., 'jungfrau', 'epix10k2M', etc.)
        out_file : str
            Path to the output CrystFEL .geom file
        """
        top = self.geo.get_top_geo()
        child = top.get_list_of_children()[0]
        x, y, z = self.geo.get_pixel_coords(oname=child.oname, oindex=0, do_tilt=True, cframe=0)
        seg = self.geo.get_seg_geo().algo
        shape = self.geo.shape3d()
        arows, acols = seg.asic_rows_cols()
        srows, _ = seg.shape()
        pix_size = seg.pixel_scale_size()
        _, nasics_in_cols = seg.number_of_asics_in_rows_cols()
        nasicsf = nasics_in_cols
        x.shape = shape
        y.shape = shape
        z.shape = shape
        distance = np.abs(np.round(np.mean(z)))*1e-6
        pyFAI_detname = psana_to_pyFAI_detname.get(detname.lower())
        txt = header_crystfel(distance=distance)
        for n in range(shape[0]):
            txt += '\n'
            for a,(r0,c0) in enumerate(seg.asic0indices()):
                vfs = np.array((\
                    x[n,r0,c0+acols-1] - x[n,r0,c0],\
                    y[n,r0,c0+acols-1] - y[n,r0,c0],\
                    z[n,r0,c0+acols-1] - z[n,r0,c0]))
                vss = np.array((\
                    x[n,r0+arows-1,c0] - x[n,r0,c0],\
                    y[n,r0+arows-1,c0] - y[n,r0,c0],\
                    z[n,r0+arows-1,c0] - z[n,r0,c0]))
                nfs = vfs/np.linalg.norm(vfs)
                nss = vss/np.linalg.norm(vss)

                pref = '\np%da%d'%(n,a)

                txt += '%s/fs = %+.6fx %+.6fy %+.6fz' % (pref, nfs[0], nfs[1], nfs[2]) \
                    + '%s/ss = %+.6fx %+.6fy %+.6fz' % (pref, nss[0], nss[1], nss[2]) \
                    + '%s/res = %.3f' % (pref, 1e6/pix_size) \
                    + '%s/corner_x = %.6f' % (pref, x[n,r0,c0]/pix_size) \
                    + '%s/corner_y = %.6f' % (pref, y[n,r0,c0]/pix_size) \
                    + '%s/coffset = %.6f' % (pref, z[n,r0,c0]*1e-6) \
                    + '%s/min_fs = %d' % (pref, (a%nasicsf)*acols) \
                    + '%s/max_fs = %d' % (pref, (a%nasicsf+1)*acols-1) \
                    + '%s/min_ss = %d' % (pref, n*srows + (a//nasicsf)*arows) \
                    + '%s/max_ss = %d' % (pref, n*srows + (a//nasicsf+1)*arows - 1) \
                    + '%s/no_index = 0' % (pref) \
                    + '\n'
                
                txt += tailer_crystfel(detname=pyFAI_detname)
        if out_file is not None:
            f = open(out_file,'w')
            f.write(txt)
            f.close()

class PsanaToPyFAI:
    """
    Class to convert psana .data geometry files directly to PyFAI detector instances

    Parameters
    ----------
    in_file : str
        Path to the input psana .data geometry file
    """
    def __init__(self, in_file, detname):
        self.geo = GeometryAccess(path=in_file, pbits=0, use_wide_pix_center=False)
        self.detector = get_detector(detname=detname, shape=self.shape3d())
        self.setup_detector()

    @classmethod
    def convert(cls, in_file, detname, image_frame=True):
        """
        Class method to convert psana .data geometry file to PyFAI detector instance

        Parameters
        ----------
        in_file : str
            Path to the input psana .data geometry file
        detname : str
            Detector name in psana (e.g., 'jungfrau', 'epix10k2M', etc.)
        image_frame : bool
            If True, use image frame coordinates; otherwise, use psana laboratory frame coordinates

        Returns
        -------
        detector : pyFAI.Detector
            PyFAI Detector instance with geometry defined by the input .data file
        """
        converter = cls(in_file, detname)
        geometry = converter.set_image_frame(image_frame=image_frame)
        converter.get_pixel_index_map(geometry=geometry, image_frame=image_frame)
        corner_array = converter.get_pixel_corners(geometry=geometry, image_frame=image_frame)
        converter.detector.set_pixel_corners(ary=corner_array)
        print(f"Converted psana .data geometry file {in_file} to PyFAI Detector.", flush=True)
        return converter.detector

    def shape3d(self):
        """
        Return 3d shape of the arrays as (<number-of-segments>, <rows>, <cols>)
        """
        seg = self.geo.get_seg_geo().algo
        sshape = seg.shape()
        ssize = seg.size()
        x,_,_ = self.geo.get_pixel_coords()
        dsize = x.size
        return (int(dsize/ssize), sshape[0], sshape[1])

    def setup_detector(self):
        """
        Pass the detector segmentation and geometry info to the PyFAI detector instance
        """
        self.detector.geo = self.geo
        self.detector.seg_geo = self.geo.get_seg_geo()
        self.detector.segname = self.geo.get_seg_geo().oname

    def set_image_frame(self, image_frame=True):
        """
        Modify the geometry to be in the image frame for PyFAI calibration

        Parameters
        ----------
        image_frame : bool
            If True, use image frame coordinates; otherwise, use psana laboratory frame coordinates
        """
        # Set reference frame to be the image frame (i.e. no offsets and undo 90° rotation)
        # X-axis: horizontal from left to right
        # Y-axis: vertical from bottom to top
        # Z-axis: opposite beam direction
        top = self.geo.get_top_geo()
        geometry = top.get_list_of_children()[0]
        self.detector.distance = -geometry.z0 * 1e-6
        if image_frame:
            geometry.x0 = 0
            geometry.y0 = 0
            geometry.z0 = 0
            geometry.rot_x = 0
            geometry.rot_y = 0
            geometry.rot_z -= 90
        return geometry

    def get_pixel_index_map(self, geometry, image_frame=True):
        """
        Create a pixel index map for assembling the detector

        Parameters
        ----------
        geometry : psana geometry object
            The geometry object from which to extract pixel coordinates
        image_frame : bool
            If True, x, y, z are image frame coordinates; otherwise, x, y, z are psana laboratory frame coordinates
        """
        x, y, _ = geometry.get_pixel_coords(do_tilt=True)
        pix_size = self.detector.pixel_size_um
        if image_frame:
            # Pixel (0,0) is at (xmin, ymax) in image coordinates
            xmin, ymax = x.min(), y.max()
            xmin, ymax = xmin - pix_size/2, ymax + pix_size/2
            rows, cols = np.array((ymax - y)/pix_size, dtype=np.uint), np.array((x - xmin)/pix_size, dtype=np.uint)
        else:
            # Pixel (0,0) is at (ymin, xmin) in psana-lab coordinates
            xmin, ymin = x.min(), y.min()
            xmin, ymin = xmin - pix_size/2, ymin - pix_size/2
            rows, cols = np.array((x - xmin)/pix_size, dtype=np.uint), np.array((y - ymin)/pix_size, dtype=np.uint)
        pixel_index_map = np.zeros((np.array(rows).shape) + (2,))
        pixel_index_map[..., 0] = rows
        pixel_index_map[..., 1] = cols
        self.detector.pixel_index_map = pixel_index_map.astype(np.int64)

    def get_pixel_corners(self, geometry, image_frame=True):
        """
        Compute the pixel corner coordinates to instantiate a 3D PyFAI detector

        Parameters
        ----------
        geometry : psana geometry object
            The geometry object from which to extract pixel coordinates
        image_frame : bool
            If True, x, y, z are image frame coordinates; otherwise, x, y, z are psana laboratory frame coordinates
        """
        x, y, z = geometry.get_pixel_coords(do_tilt=True)
        x, y, z = image_to_pyfai(x, y, z, image_frame)
        npanels = self.detector.n_modules
        nasics = self.detector.n_asics
        asics_shape = self.detector.asics_shape
        fs_size = self.detector.fs_size
        ss_size = self.detector.ss_size
        pixel_corners = np.zeros([npanels * ss_size * asics_shape[0], fs_size * asics_shape[1], 4, 3])
        for p in range(npanels):
            xp = x[p, :]
            yp = y[p, :]
            zp = z[p, :]
            vfs = np.array((\
                xp[0, fs_size * asics_shape[1] - 1] - xp[0, 0],\
                yp[0, fs_size * asics_shape[1] - 1] - yp[0, 0],\
                zp[0, fs_size * asics_shape[1] - 1] - zp[0, 0]))
            vss = np.array((\
                xp[ss_size * asics_shape[0] - 1,0] - xp[0, 0],\
                yp[ss_size * asics_shape[0] - 1,0] - yp[0, 0],\
                zp[ss_size * asics_shape[0] - 1,0] - zp[0, 0]))
            nfs = vfs / np.linalg.norm(vfs)
            nss = vss / np.linalg.norm(vss)
            for a in range(nasics):
                if nasics == 1:
                    arow = 0
                    acol = 0
                else:
                    arow = a // (nasics//2)
                    acol = a % (nasics//2)
                fs_portion = slice(acol * fs_size, (acol + 1) * fs_size)
                ss_portion = slice(arow * ss_size, (arow + 1) * ss_size)
                slab_offset = p * asics_shape[0] * ss_size
                fs_portion_slab = slice(acol * fs_size, (acol + 1) * fs_size)
                ss_portion_slab = slice(arow * ss_size + slab_offset, (arow + 1) * ss_size + slab_offset)
                ssx, ssy, ssz = nss * self.detector.pixel1
                fsx, fsy, fsz = nfs * self.detector.pixel2
                xasic = x[p, ss_portion, fs_portion]
                yasic = y[p, ss_portion, fs_portion]
                zasic = z[p, ss_portion, fs_portion]
                # psana reports pixel centers, so the corners sit half a pixel either side of
                # them. Offsets of 0 and 1 would put corner 0 on the center and shift the whole
                # grid outward by half a pixel.
                corner_ss = np.array([-0.5, -0.5, 0.5, 0.5])
                corner_fs = np.array([-0.5, 0.5, 0.5, -0.5])
                xasic = xasic[:, :, np.newaxis] + corner_ss * ssx + corner_fs * fsx
                yasic = yasic[:, :, np.newaxis] + corner_ss * ssy + corner_fs * fsy
                zasic = zasic[:, :, np.newaxis] + corner_ss * ssz + corner_fs * fsz
                pixel_corners[ss_portion_slab, fs_portion_slab, :, 0] = zasic
                pixel_corners[ss_portion_slab, fs_portion_slab, :, 1] = xasic
                pixel_corners[ss_portion_slab, fs_portion_slab, :, 2] = yasic
        return pixel_corners

class PsanaToDXTBX:
    """
    Class to convert psana .data geometry files to dxtbx Detector models

    Emits one flat panel per ASIC, with absolute lab-frame origins in millimeters. This is the
    form consumed by simtbx.nanoBragg: panels map 1:1 onto the psana calibration array, so no
    pixel padding is required, and inter-ASIC gaps are carried by the panel origins.

    All detector placement is zeroed (the camera node offsets, its 90 degree rotation, and the
    sample-detector distance), so the panel origins hold pure metrology and every origin has
    z = 0. Placement is expected to come from geometric parameters applied downstream. Because
    the total distance is therefore zero, ``get_beam_centre()`` and ``get_ray_intersection()``
    raise on the returned model; use PsanaToCCTBX if an absolutely placed detector is needed.
    The distance read from the input file is still available as ``self.detz_mm``.

    Parameters
    ----------
    in_file : str
        Path to the input psana .data geometry file
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
    """
    def __init__(self, in_file, detname):
        self.geo = GeometryAccess(path=in_file, pbits=0, use_wide_pix_center=False)
        check_detname(detname=detname, shape=self.shape3d())
        self.detz_mm = None

    @classmethod
    def convert(cls, in_file, detname):
        """
        Class method to convert psana .data geometry file to a dxtbx Detector

        Parameters
        ----------
        in_file : str
            Path to the input psana .data geometry file
        detname : str
            Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)

        Returns
        -------
        detector : dxtbx.model.Detector
            dxtbx Detector with one flat panel per ASIC, origins in millimeters
        """
        from dxtbx.model import DetectorFactory

        converter = cls(in_file, detname)
        panels = converter.get_panels()
        detector = DetectorFactory.from_dict({"panels": panels})
        print(f"Converted psana .data geometry file {in_file} to dxtbx Detector.", flush=True)
        return detector

    @classmethod
    def convert_to_file(cls, in_file, detname, out_file):
        """
        Class method to convert psana .data geometry file to a dxtbx .expt geometry file

        Parameters
        ----------
        in_file : str
            Path to the input psana .data geometry file
        detname : str
            Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
        out_file : str
            Path to the output dxtbx .expt/.json geometry file
        """
        detector = cls.convert(in_file=in_file, detname=detname)
        write_experiment_list(detector=detector, out_file=out_file)
        return detector

    def shape3d(self):
        """
        Return 3d shape of the arrays as (<number-of-segments>, <rows>, <cols>)
        """
        seg = self.geo.get_seg_geo().algo
        sshape = seg.shape()
        ssize = seg.size()
        x,_,_ = self.geo.get_pixel_coords()
        dsize = x.size
        return (int(dsize/ssize), sshape[0], sshape[1])

    def zero_placement(self):
        """
        Zero the detector placement so that pixel coordinates hold metrology only

        Removes the camera node offsets and the 90 degree rotation of the psana convention,
        and records the sample-detector distance in millimeters on ``self.detz_mm``.

        Returns
        -------
        geometry : GeometryObject
            The camera node, with its placement zeroed
        """
        top = self.geo.get_top_geo()
        geometry = top.get_list_of_children()[0]
        self.detz_mm = -geometry.z0 * 1e-3
        geometry.x0 = 0
        geometry.y0 = 0
        geometry.z0 = 0
        geometry.rot_x = 0
        geometry.rot_y = 0
        geometry.rot_z -= 90
        return geometry

    def get_panels(self):
        """
        Build one dxtbx panel dictionary per ASIC

        psana reports pixel centers in microns; dxtbx expects the outer corner of the first
        pixel in millimeters, so coordinates are scaled by 1e-3 and shifted by half a pixel
        along both basis vectors.

        Returns
        -------
        panels : list of dict
            Panel dictionaries suitable for ``DetectorFactory.from_dict``
        """
        geometry = self.zero_placement()
        shape = self.shape3d()
        seg = self.geo.get_seg_geo().algo
        arows, acols = seg.asic_rows_cols()
        pix = seg.pixel_scale_size() * 1e-3
        x, y, z = geometry.get_pixel_coords(do_tilt=True)
        x.shape = y.shape = z.shape = shape
        panels = []
        for n in range(shape[0]):
            for a, (r0, c0) in enumerate(seg.asic0indices()):
                vfs = np.array((\
                    x[n,r0,c0+acols-1] - x[n,r0,c0],\
                    y[n,r0,c0+acols-1] - y[n,r0,c0],\
                    z[n,r0,c0+acols-1] - z[n,r0,c0]))
                vss = np.array((\
                    x[n,r0+arows-1,c0] - x[n,r0,c0],\
                    y[n,r0+arows-1,c0] - y[n,r0,c0],\
                    z[n,r0+arows-1,c0] - z[n,r0,c0]))
                nfs = vfs / np.linalg.norm(vfs)
                nss = vss / np.linalg.norm(vss)
                center = np.array((x[n,r0,c0], y[n,r0,c0], z[n,r0,c0])) * 1e-3
                origin = center - 0.5 * pix * nfs - 0.5 * pix * nss
                panels.append({
                    "name": "ARRAY_D0M%dA%d" % (n, a),
                    "type": "SENSOR_PAD",
                    "fast_axis": tuple(nfs),
                    "slow_axis": tuple(nss),
                    "origin": tuple(origin),
                    "raw_image_offset": (0, 0),
                    "image_size": (acols, arows),
                    "pixel_size": (pix, pix),
                    "trusted_range": (-10.0, 2e6),
                    "thickness": 0.0,
                    "material": "",
                    "mu": 0.0,
                    "identifier": "",
                    "mask": [],
                    "gain": 1.0,
                    "pedestal": 0.0,
                    "px_mm_strategy": {"type": "SimplePxMmStrategy"},
                })
        return panels

class PsanaToCCTBX:
    """
    Class to convert psana .data geometry files to cctbx-convention dxtbx Detector models

    Emits the geometry the way cctbx.xfel does, reproducing the structure of refined .expt
    files: one padded module-level panel per psana segment, arranged in a three level
    hierarchy (detector -> module groups -> panels) whose local frames compose to the absolute
    positions. Nothing is zeroed, so the sample-detector distance from the input file is kept.

    Panels are padded to the cctbx pixel layout, which inserts the double-width pixels at the
    inter-ASIC boundaries (a 512x1024 Jungfrau segment becomes a 514x1030 panel). Raw psana
    frames must therefore be passed through ``serialtbx.detector.jungfrau.correct_panel``
    before they can be indexed with this model.

    Parameters
    ----------
    in_file : str
        Path to the input psana .data geometry file
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
    """
    def __init__(self, in_file, detname):
        self.geo = GeometryAccess(path=in_file, pbits=0, use_wide_pix_center=False)
        check_detname(detname=detname, shape=self.shape3d())

    @classmethod
    def convert(cls, in_file, detname):
        """
        Class method to convert psana .data geometry file to a cctbx-convention dxtbx Detector

        Parameters
        ----------
        in_file : str
            Path to the input psana .data geometry file
        detname : str
            Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)

        Returns
        -------
        detector : dxtbx.model.Detector
            dxtbx Detector with a cctbx-style hierarchy and padded module-level panels
        """
        converter = cls(in_file, detname)
        detector = converter.build_detector()
        print(f"Converted psana .data geometry file {in_file} to cctbx dxtbx Detector.", flush=True)
        return detector

    @classmethod
    def convert_to_file(cls, in_file, detname, out_file):
        """
        Class method to convert psana .data geometry file to a cctbx-convention geometry file

        Parameters
        ----------
        in_file : str
            Path to the input psana .data geometry file
        detname : str
            Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
        out_file : str
            Path to the output dxtbx .expt/.json geometry file
        """
        detector = cls.convert(in_file=in_file, detname=detname)
        write_experiment_list(detector=detector, out_file=out_file)
        return detector

    def shape3d(self):
        """
        Return 3d shape of the arrays as (<number-of-segments>, <rows>, <cols>)
        """
        seg = self.geo.get_seg_geo().algo
        sshape = seg.shape()
        ssize = seg.size()
        x,_,_ = self.geo.get_pixel_coords()
        dsize = x.size
        return (int(dsize/ssize), sshape[0], sshape[1])

    def build_detector(self):
        """
        Build the dxtbx Detector hierarchy from the psana geometry tree

        Returns
        -------
        detector : dxtbx.model.Detector
            dxtbx Detector with a cctbx-style hierarchy
        """
        from dxtbx.model import Detector
        from scitbx.matrix import col
        from serialtbx.detector.xtc import basis_from_geo
        try:
            from psana.pscalib.geometry.SegGeometryStore import sgs  # psana2
        except ImportError:
            from PSCalib.SegGeometryStore import sgs  # psana1

        pix = self.geo.get_pixel_scale_size() * 1e-3
        detector = Detector()
        group = detector.hierarchy()

        # Collapse the single-child chain at the top of the tree (e.g. IP -> CAMERA)
        root = self.geo.get_top_geo()
        root_basis = basis_from_geo(root)
        while len(root.get_list_of_children()) == 1:
            child = root.get_list_of_children()[0]
            root_basis = root_basis * basis_from_geo(child)
            root = child

        # The detector sits downstream of the sample, so its origin has negative z
        translation = root_basis.translation
        root_basis.translation = col((translation[0], translation[1], -abs(translation[2])))
        self.set_local_frame(node=group, basis=root_basis)
        group.set_name("D0")

        for n, module in enumerate(root.get_list_of_children()):
            subgroup = group.add_group()
            self.set_local_frame(node=subgroup, basis=basis_from_geo(module))
            subgroup.set_name("D0M%d" % n)

            seg_geo = sgs.Create(segname=module.oname)
            xx, yy = seg_geo.get_seg_xy_maps_um()
            xx = xx * 1e-3
            yy = yy * 1e-3
            arows, acols = self.panel_rows_cols(seg_geo)
            for a, (r0, c0) in enumerate(self.panel0indices(seg_geo)):
                # cctbx takes the first pixel center directly as the panel origin
                origin = np.array((xx[r0][c0], yy[r0][c0], 0.0))
                vfs = np.array((xx[r0][c0+1], yy[r0][c0+1], 0.0)) - origin
                vss = np.array((xx[r0+1][c0], yy[r0+1][c0], 0.0)) - origin
                panel = subgroup.add_panel()
                panel.set_local_frame(
                    tuple(vfs / np.linalg.norm(vfs)),
                    tuple(vss / np.linalg.norm(vss)),
                    tuple(origin),
                )
                panel.set_name("ARRAY_D0M%dA%d" % (n, a))
                panel.set_type("SENSOR_PAD")
                panel.set_pixel_size((pix, pix))
                panel.set_image_size((acols, arows))
                panel.set_trusted_range((-10.0, 2e6))
        return detector

    def is_padded(self, seg_geo):
        """
        Whether the segment spans a whole number of pixels once its wide pixels are counted

        Jungfrau segments have double-width pixels at the inter-ASIC boundaries that span an
        integral number of pixel pitches, so cctbx represents the whole segment as one padded
        panel (512x1024 -> 514x1030). Other detectors, such as ePix, have fractional gaps that
        cannot be padded this way, and cctbx emits one panel per ASIC instead.

        Parameters
        ----------
        seg_geo : SegGeometry
            Segment geometry object from the segment geometry store

        Returns
        -------
        padded : bool
            True if the segment can be represented as a single padded panel
        """
        rows, cols = self.padded_rows_cols(seg_geo)
        return abs(rows - round(rows)) < 1e-6 and abs(cols - round(cols)) < 1e-6

    def padded_rows_cols(self, seg_geo):
        """
        Extent of the segment in pixel pitches, counting its wide pixels

        Parameters
        ----------
        seg_geo : SegGeometry
            Segment geometry object from the segment geometry store

        Returns
        -------
        rows, cols : float, float
            Slow and fast extent of the segment, in units of the pixel pitch
        """
        xx, yy = seg_geo.get_seg_xy_maps_um()
        pix = seg_geo.pixel_scale_size()
        cols = abs(xx[0, -1] - xx[0, 0]) / pix + 1
        rows = abs(yy[-1, 0] - yy[0, 0]) / pix + 1
        return rows, cols

    def panel_rows_cols(self, seg_geo):
        """
        Shape of a single cctbx panel for this segment

        Returns
        -------
        rows, cols : int, int
            Slow and fast size of the panel, in pixels
        """
        if self.is_padded(seg_geo):
            rows, cols = self.padded_rows_cols(seg_geo)
            return int(round(rows)), int(round(cols))
        return seg_geo.asic_rows_cols()

    def panel0indices(self, seg_geo):
        """
        Index of the first pixel of each cctbx panel within the segment

        Returns
        -------
        indices : tuple
            One (row, col) pair per panel emitted for this segment
        """
        if self.is_padded(seg_geo):
            return ((0, 0),)
        return seg_geo.asic0indices()

    def set_local_frame(self, node, basis):
        """
        Apply a cctbx basis to a dxtbx hierarchy node as its local frame

        Parameters
        ----------
        node : dxtbx.model.Panel or dxtbx.model.DetectorNode
            Node whose local frame is set
        basis : serialtbx.detector.basis
            Basis holding the orientation and translation of the node
        """
        from scitbx.matrix import col

        origin = col((basis * col((0, 0, 0, 1)))[0:3])
        fast = col((basis * col((1, 0, 0, 1)))[0:3]) - origin
        slow = col((basis * col((0, 1, 0, 1)))[0:3]) - origin
        node.set_local_frame(fast.elems, slow.elems, origin.elems)

class PyFAIToPsana:
    """
    Class to convert PyFAI .poni files to psana .data geometry files

    Parameters
    ----------
    in_file : str
        Path to the PyFAI .poni file containing detector geometry parameters
    detector : pyFAI.Detector
        PyFAI Detector model to be calibrated by the .poni file
    """
    def __init__(self, in_file, detector):
        self.detector = detector
        ai = pyFAI.load(in_file)
        self.params = ai.param
        self.X, self.Y, self.Z = correct_geometry(detector=self.detector, params=self.params)

    @classmethod
    def convert(cls, in_file, detector, out_file, image_frame=True):
        """
        Class method to convert PyFAI .poni file to psana .data geometry file
        
        Parameters
        ----------
        in_file : str
            Path to the PyFAI .poni file containing detector geometry parameters
        detector : pyFAI.Detector
            PyFAI Detector model to be calibrated by the .poni file
        out_file : str
            Path to the output psana .data geometry file
        image_frame : bool
            If True, use image frame coordinates; otherwise, use psana laboratory frame coordinates
        """
        converter = cls(in_file, detector)
        converter.convert_pyfai_to_data(out_file=out_file, image_frame=image_frame)
        print(f"Converted PyFAI .poni file {in_file} to psana .data geometry file {out_file}.", flush=True)

    def convert_pyfai_to_data(self, out_file, image_frame=True):
        """
        Main function to convert PyFAI coordinates to psana .data geometry file

        Parameters
        ----------
        out_file : str
            Path to the output .data file
        image_frame : bool
            If True, use image frame coordinates; otherwise, use psana laboratory frame coordinates
        """
        x = self.X.reshape(self.detector.calib_shape)
        y = self.Y.reshape(self.detector.calib_shape)
        z = self.Z.reshape(self.detector.calib_shape)
        x, y, z = pyfai_to_image(x, y, z)
        geo = self.detector.geo
        top = geo.get_top_geo()
        child = top.get_list_of_children()[0]
        topname = top.oname
        childname = child.oname
        npanels = self.detector.n_modules
        asics_shape = self.detector.asics_shape
        fs_size = self.detector.fs_size
        ss_size = self.detector.ss_size
        recs = header_psana(detname=self.detector.detname)
        distance = round(np.mean(z))
        for p in range(npanels):
            xp = x[p, :]
            yp = y[p, :]
            zp = z[p, :]
            vfs = np.array((\
                xp[0, fs_size * asics_shape[1] - 1] - xp[0, 0],\
                yp[0, fs_size * asics_shape[1] - 1] - yp[0, 0],\
                zp[0, fs_size * asics_shape[1] - 1] - zp[0, 0]))
            vss = np.array((\
                xp[ss_size * asics_shape[0] - 1,0] - xp[0, 0],\
                yp[ss_size * asics_shape[0] - 1,0] - yp[0, 0],\
                zp[ss_size * asics_shape[0] - 1,0] - zp[0, 0]))
            nfs = vfs / np.linalg.norm(vfs)
            nss = vss / np.linalg.norm(vss)
            vcent = (np.mean(xp), np.mean(yp), np.mean(zp)-distance)
            angle_deg_z = degrees(atan2(nfs[1], nfs[0]))
            angle_z, tilt_z = angle_and_tilt(angle_deg_z)
            tilt_x, tilt_y = tilt_xy(nfs, nss)
            angle_x = child.get_list_of_children()[p].rot_x
            angle_y = child.get_list_of_children()[p].rot_y
            recs += '\n%12s  0 %12s %2d' %(childname, self.detector.segname, p)\
                +'  %8d %8d %8d %7.0f %6.0f %6.0f   %8.5f  %8.5f  %8.5f'%\
                (vcent[0], vcent[1], vcent[2], angle_z, angle_y, angle_x, tilt_z, tilt_y, tilt_x)
        if image_frame:
            recs += '\n%12s  0 %12s  0' %(topname, childname)\
                +'         0        0 %8d       0      0      0    0.00000   0.00000   0.00000' % (distance)
        else:
            recs += '\n%12s  0 %12s  0' %(topname, childname)\
                +'         0        0 %8d      90      0      0    0.00000   0.00000   0.00000' % (distance)
        f=open(out_file,'w')
        f.write(recs)
        f.close()

class PyFAIToCrystFEL:
    """
    Class to write CrystFEL .geom geometry files from PyFAI .poni files

    Parameters
    ----------
    in_file : str
        Path to the input .poni file containing detector geometry parameters
    detector : pyFAI.Detector
        PyFAI detector model to be calibrated by the .poni file
    """

    def __init__(self, in_file, detector):
        self.detector = detector
        ai = pyFAI.load(in_file)
        self.params = ai.param
        self.X, self.Y, self.Z = correct_geometry(detector=self.detector, params=self.params)

    @classmethod
    def convert(cls, in_file, detector, out_file, image_frame=True):
        """
        Class method to convert PyFAI .poni file to CrystFEL .geom geometry file

        Parameters
        ----------
        in_file : str
            Path to the input .poni file containing detector geometry parameters
        detector : pyFAI.Detector
            PyFAI detector model to be calibrated by the .poni file
        out_file : str
            Path to the output CrystFEL .geom file
        """
        converter = cls(in_file, detector)
        converter.convert_pyfai_to_geom(out_file=out_file, image_frame=image_frame)
        print(f"Converted PyFAI .poni file {in_file} to CrystFEL .geom geometry file {out_file}.")
    
    def convert_pyfai_to_geom(self, out_file, image_frame=True):
        """
        Main function to convert PyFAI coordinates to CrystFEL .geom geometry file

        Parameters
        ----------
        output_file : str
            Path to the output .geom file
        """
        x = self.X.reshape(self.detector.calib_shape)
        y = self.Y.reshape(self.detector.calib_shape)
        z = self.Z.reshape(self.detector.calib_shape)
        x, y, z = pyfai_to_image(x, y, z, image_frame)
        seg = self.detector.seg_geo.algo
        nsegs = int(x.size/seg.size())
        arows, acols = seg.asic_rows_cols()
        srows, _ = seg.shape()
        pix_size = seg.pixel_scale_size()
        _, nasics_in_cols = seg.number_of_asics_in_rows_cols()
        nasicsf = nasics_in_cols
        txt = header_crystfel(self.params[0])
        for n in range(nsegs):
            txt += '\n'
            for a,(r0,c0) in enumerate(seg.asic0indices()):
                vfs = np.array((\
                    x[n,r0,c0+acols-1] - x[n,r0,c0],\
                    y[n,r0,c0+acols-1] - y[n,r0,c0],\
                    z[n,r0,c0+acols-1] - z[n,r0,c0]))
                vss = np.array((\
                    x[n,r0+arows-1,c0] - x[n,r0,c0],\
                    y[n,r0+arows-1,c0] - y[n,r0,c0],\
                    z[n,r0+arows-1,c0] - z[n,r0,c0]))
                nfs = vfs/np.linalg.norm(vfs)
                nss = vss/np.linalg.norm(vss)

                pref = '\np%da%d'%(n,a)

                txt +='%s/fs = %+.6fx %+.6fy %+.6fz' % (pref, nfs[0], nfs[1], nfs[2])\
                    + '%s/ss = %+.6fx %+.6fy %+.6fz' % (pref, nss[0], nss[1], nss[2])\
                    + '%s/res = %.3f' % (pref, 1e6/pix_size)\
                    + '%s/corner_x = %.6f' % (pref, x[n,r0,c0]/pix_size)\
                    + '%s/corner_y = %.6f' % (pref, y[n,r0,c0]/pix_size)\
                    + '%s/coffset = %.6f' % (pref, z[n,r0,c0]*1e-6)\
                    + '%s/min_fs = %d' % (pref, (a%nasicsf)*acols)\
                    + '%s/max_fs = %d' % (pref, (a%nasicsf+1)*acols-1)\
                    + '%s/min_ss = %d' % (pref, n*srows + (a//nasicsf)*arows)\
                    + '%s/max_ss = %d' % (pref, n*srows + (a//nasicsf+1)*arows - 1)\
                    + '%s/no_index = 0' % (pref)\
                    + '\n'
            
                txt += tailer_crystfel(self.detector.detname)

        if out_file is not None:
            f = open(out_file,'w')
            f.write(txt)
            f.close()

class CrystFELToPsana:
    """
    Class to convert CrystFEL .geom geometry files to psana .data geometry files

    Parameters
    ----------
    in_file : str
        Path to the CrystFEL .geom file
    detname : str
        Detector name
    """
    def __init__(self, in_file, detname):
        self.valid = False
        self.load_geom(in_file=in_file)
        template = get_default_geometry(detname=detname)
        self.geo = GeometryAccess(path=template, pbits=0, use_wide_pix_center=False)

    @classmethod
    def convert(cls, in_file, detname, out_file):
        """
        Class method to convert CrystFEL .geom geometry file to psana .data geometry file

        Parameters
        ----------
        in_file : str
            Path to the CrystFEL .geom file
        detname : str
            Detector name
        out_file : str
            Path to the output psana .data geometry file
        """
        converter = cls(in_file, detname)
        converter.convert_geom_to_data(detname=detname, out_file=out_file)
        print(f"Converted CrystFEL .geom file {in_file} to psana .data geometry file {out_file}.", flush=True)

    def _parse_line_as_parameter(self, line):
        assert isinstance(line, str), 'line is not a str object'
        fields = line.split()
        nfields = len(fields)
        if fields[1] != '=':
            self.list_of_ignored_records.append(line)
            return
        keys = fields[0].split('/')
        nkeys = len(keys)
        if nkeys==1:
            if nfields>3:
                self.list_of_ignored_records.append(line)
                return
            k0 = keys[0]
            self.dict_of_pars[k0] = float(fields[2]) if k0 in ('res', 'adu_per_eV', 'coffset') else\
                ' '.join(fields[2:])
        elif nkeys==2:
            k0, k1 = keys
            resp = str_is_segment_and_asic(k0)
            if resp: k0=resp
            v = '' if nfields<3 else\
                sfields_to_xyz_vector(fields[2:]) if k1 in ('fs','ss') else\
                int(fields[2]) if k1 in ('max_ss', 'min_ss', 'max_fs', 'min_fs', 'no_index') else\
                int(fields[2]) if k1 in ('max_x', 'min_x', 'max_y', 'min_y') else\
                float(fields[2]) if k1 in ('res', 'corner_x', 'corner_y', 'adu_per_eV', 'coffset') else\
                float(fields[2]) if k1 in ('xfs', 'yfs', 'xss', 'yss') else\
                ' '.join(fields[2:])
            if k0 in self.dict_of_pars.keys():
                self.dict_of_pars[k0][k1] = v
            else:
                self.dict_of_pars[k0] = {k1:v,}
        else:
            self.list_of_ignored_records.append(line)
            return

    def load_geom(self, in_file):
        self.valid = False
        self.list_of_comments = []
        self.list_of_ignored_records = []
        self.dict_of_pars = {}
        f=open(in_file,'r')
        for linef in f:
            line = linef.strip('\n')
            if not line.strip(): continue
            if line[0] == ';':
                self.list_of_comments.append(line)
                continue
            self._parse_line_as_parameter(line)
        f.close()
        self.valid = True

    def geom_to_data(self, panelasics, detname, out_file):
        geo = self.geo
        top = geo.get_top_geo()
        child = top.get_list_of_children()[0]
        topname = top.oname
        childname = child.oname
        seg = geo.get_seg_geo()
        sg = seg.algo
        segname = seg.oname
        X, Y, Z = sg.pixel_coord_array()
        PIX_SIZE_UM = sg.get_pix_size_um()
        M_TO_UM = 1e6
        xc0, yc0, _ = X[0,0], Y[0,0], Z[0,0]
        distance = self.dict_of_pars.get('coffset', 0)
        recs = header_psana(detname=detname)
        segz = np.array([self.dict_of_pars[k].get('coffset', 0) for k in panelasics.split(',')])
        meanroundz = round(segz.mean()*1e6)*1e-6
        distance += meanroundz
        for p, panel in enumerate(panelasics.split(',')):
            dicasic = self.dict_of_pars[panel]
            nfs = np.array(dicasic.get('fs', None), dtype=np.float64) 
            nss = np.array(dicasic.get('ss', None), dtype=np.float64) 
            vfs = nfs*abs(xc0)
            vss = nfs*abs(yc0)
            x0pix = dicasic.get('corner_x', 0)
            y0pix = dicasic.get('corner_y', 0)
            z0    = dicasic.get('coffset', 0)
            v00center = vfs + vss
            v00corner = np.array((x0pix*PIX_SIZE_UM, y0pix*PIX_SIZE_UM, (z0 - distance)*M_TO_UM))
            vcent = v00corner + v00center
            angle_deg_z = degrees(atan2(nfs[1], nfs[0]))
            angle_z, tilt_z = angle_and_tilt(angle_deg_z)
            tilt_x, tilt_y = tilt_xy(nfs, nss)
            angle_x = child.get_list_of_children()[p].rot_x
            angle_y = child.get_list_of_children()[p].rot_y
            recs += '\n%12s  0 %12s %2d' %(childname, segname, p)\
                +'  %8d %8d %8d %7.0f %6.0f %6.0f   %8.5f  %8.5f  %8.5f'%\
                (vcent[0], vcent[1], vcent[2], angle_z, angle_y, angle_x, tilt_z, tilt_y, tilt_x)
        recs += '\n%12s  0 %12s  0' %(topname, childname)\
            +'         0        0 %8d       0      0      0    0.00000   0.00000   0.00000' % (distance*M_TO_UM)
        f=open(out_file,'w')
        f.write(recs)
        f.close()

    def convert_geom_to_data(self, detname, out_file):
        pyFAI_detname = psana_to_pyFAI_detname.get(detname.lower())
        if pyFAI_detname is None:
            raise ValueError(f"Detector name {detname} not recognized.")
        panelasics = detname_to_pars.get(pyFAI_detname.lower())
        if panelasics is None:
            raise ValueError(f"Detector name {detname} not implemented.")
        self.geom_to_data(panelasics, pyFAI_detname, out_file)

class CrystFELToPyFAI:
    """
    Class to convert CrystFEL .geom geometry files to PyFAI Detector model by using intermediate psana .data files

    Parameters
    ----------
    in_file : str
        Path to the CrystFEL .geom file
    """
    def __init__(self, in_file):
        path = os.path.dirname(in_file)
        self.temp_file = os.path.join(path, "temp.data")
    
    @classmethod
    def convert(cls, in_file, detname):
        """
        Class method to convert CrystFEL .geom geometry file to PyFAI Detector instance

        Parameters
        ----------
        in_file : str
            Path to the CrystFEL .geom file
        detname : str
            Detector name

        Returns
        -------
        detector : pyFAI.Detector
            PyFAI Detector instance with geometry defined by the input .geom file
        """
        converter = cls(in_file)
        CrystFELToPsana.convert(in_file=in_file, detname=detname, out_file=converter.temp_file)
        detector = PsanaToPyFAI.convert(in_file=converter.temp_file, detname=detname)
        os.unlink(converter.temp_file)
        print(f"Converted CrystFEL .geom file {in_file} to PyFAI Detector instance.", flush=True)
        return detector
"""
Converter classes to convert between psana .data geometry files, CrystFEL .geom geometry files,
PyFAI Detector models, and dxtbx Detector models
"""

try:
    from psana.pscalib.geometry.GeometryAccess import GeometryAccess  # psana2
except ImportError:
    from PSCalib.GeometryAccess import GeometryAccess  # psana1

import os
import tempfile
from pathlib import Path
import numpy as np
from math import atan2, degrees

from LCLSGeom.manager import get_default_geometry
from LCLSGeom.utils.crystfel_utils import (
    header_crystfel,
    panel_names,
    panel_crystfel,
    sfields_to_xyz_vector, 
    str_is_segment_and_asic,
    tailer_crystfel
)
from LCLSGeom.utils.detectors import check_detname
from LCLSGeom.utils.dxtbx_utils import panel_layout, set_local_frame, write_experiment
from LCLSGeom.utils.geomcalc import asic_axes
from LCLSGeom.utils.pscalib_utils import (angle_and_tilt, get_root_geometry, header_psana,
                                          shape3d, tilt_xy)
from LCLSGeom.utils.pyfai_utils import (correct_geometry, get_pyFAI_detname, image_to_pyfai,
                                        load_poni, pyfai_to_image)


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
        shape = shape3d(self.geo)
        arows, acols = seg.asic_rows_cols()
        srows, _ = seg.shape()
        pix_size = seg.pixel_scale_size()
        _, nasicsf = seg.number_of_asics_in_rows_cols()
        x.shape = shape
        y.shape = shape
        z.shape = shape
        distance = np.abs(np.round(np.mean(z)))*1e-6
        pyFAI_detname = get_pyFAI_detname(detname)
        txt = header_crystfel(distance=distance)
        for n in range(shape[0]):
            txt += '\n'
            for a,(r0,c0) in enumerate(seg.asic0indices()):
                nfs, nss = asic_axes(x[n], y[n], z[n], r0, c0, arows, acols)
                txt += panel_crystfel(n, a, nfs, nss, x[n,r0,c0], y[n,r0,c0], z[n,r0,c0],
                                      pix_size, arows, acols, srows, nasicsf)
                txt += tailer_crystfel(detname=pyFAI_detname)
        Path(out_file).write_text(txt)

class PsanaToPyFAI:
    """
    Class to convert psana .data geometry files directly to PyFAI detector instances

    Parameters
    ----------
    in_file : str
        Path to the input psana .data geometry file
    detname : str
        Detector name in psana (e.g., 'jungfrau', 'epix10k2M', etc.)
    """
    def __init__(self, in_file, detname):
        from LCLSGeom.detector import get_detector

        self.geo = GeometryAccess(path=in_file, pbits=0, use_wide_pix_center=False)
        self.detector = get_detector(detname=detname, shape=shape3d(self.geo))
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
        camera = converter.geo.get_top_geo().get_list_of_children()[0]
        converter.detector.distance = -camera.z0 * 1e-6
        geometry = get_root_geometry(geo=converter.geo, image_frame=image_frame)
        converter.get_pixel_index_map(geometry=geometry, image_frame=image_frame)
        corner_array = converter.get_pixel_corners(geometry=geometry, image_frame=image_frame)
        converter.detector.set_pixel_corners(ary=corner_array)
        print(f"Converted psana .data geometry file {in_file} to PyFAI Detector.", flush=True)
        return converter.detector

    def setup_detector(self):
        """
        Pass the detector segmentation and geometry info to the PyFAI detector instance
        """
        self.detector.geo = self.geo
        self.detector.seg_geo = self.geo.get_seg_geo()
        self.detector.segname = self.geo.get_seg_geo().oname

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
            # Pixel (0,0) is at (xmin, ymin) in psana-lab coordinates, with rows along x
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

        Returns
        -------
        pixel_corners : np.ndarray
            Corner coordinates in meters, shaped (<slab-rows>, <slab-cols>, 4, 3) with the last
            axis ordered (z, x, y) as PyFAI expects
        """
        x, y, z = geometry.get_pixel_coords(do_tilt=True)
        x, y, z = image_to_pyfai(x, y, z, image_frame)
        npanels = self.detector.n_modules
        nasics = self.detector.n_asics
        asics_shape = self.detector.asics_shape
        fs_size = self.detector.fs_size
        ss_size = self.detector.ss_size
        pixel_corners = np.zeros([npanels * ss_size * asics_shape[0], fs_size * asics_shape[1], 4, 3])
        # psana reports pixel centers, so the corners sit half a pixel either side of them.
        # Offsets of 0 and 1 would put corner 0 on the center and shift the whole grid outward
        # by half a pixel.
        corner_ss = np.array([-0.5, -0.5, 0.5, 0.5])
        corner_fs = np.array([-0.5, 0.5, 0.5, -0.5])
        for p in range(npanels):
            nfs, nss = asic_axes(x[p], y[p], z[p], 0, 0,
                                 ss_size * asics_shape[0], fs_size * asics_shape[1])
            for a in range(nasics):
                arow, acol = a // asics_shape[1], a % asics_shape[1]
                fs_portion = slice(acol * fs_size, (acol + 1) * fs_size)
                ss_portion = slice(arow * ss_size, (arow + 1) * ss_size)
                # Panels are stacked along the slow axis of the PyFAI slab
                slab_offset = p * asics_shape[0] * ss_size
                ss_portion_slab = slice(ss_portion.start + slab_offset, ss_portion.stop + slab_offset)
                ssx, ssy, ssz = nss * self.detector.pixel1
                fsx, fsy, fsz = nfs * self.detector.pixel2
                xasic = x[p, ss_portion, fs_portion]
                yasic = y[p, ss_portion, fs_portion]
                zasic = z[p, ss_portion, fs_portion]
                xasic = xasic[:, :, np.newaxis] + corner_ss * ssx + corner_fs * fsx
                yasic = yasic[:, :, np.newaxis] + corner_ss * ssy + corner_fs * fsy
                zasic = zasic[:, :, np.newaxis] + corner_ss * ssz + corner_fs * fsz
                pixel_corners[ss_portion_slab, fs_portion, :,0] = zasic
                pixel_corners[ss_portion_slab, fs_portion, :,1] = xasic
                pixel_corners[ss_portion_slab, fs_portion, :,2] = yasic
        return pixel_corners

class PsanaToDXTBX:
    """
    Class to convert psana .data geometry files to flat dxtbx Detector models.

    Use PsanaToCCTBX instead when an absolutely placed detector is needed, or when the consumer
    expects the cctbx.xfel panel layout and hierarchy.

    ``convert`` returns the Detector model.
    ``convert_to_file`` writes the same geometry as an .expt file.

    Parameters
    ----------
    in_file : str
        Path to the input psana .data geometry file
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
    """
    def __init__(self, in_file, detname):
        self.geo = GeometryAccess(path=in_file, pbits=0, use_wide_pix_center=False)
        check_detname(detname=detname, shape=shape3d(self.geo))

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
        detector = cls(in_file, detname).build_detector()
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

        Returns
        -------
        detector : dxtbx.model.Detector
            The Detector that was written
        """
        detector = cls.convert(in_file=in_file, detname=detname)
        write_experiment(detector=detector, out_file=out_file)
        return detector

    def build_detector(self):
        """
        Build a flat dxtbx Detector with one panel per ASIC

        psana reports pixel centers in microns; dxtbx expects the outer corner of the first
        pixel in millimeters, so coordinates are scaled by 1e-3 and shifted by half a pixel
        along both basis vectors.

        Returns
        -------
        detector : dxtbx.model.Detector
            dxtbx Detector with one flat panel per ASIC, origins in millimeters
        """
        from dxtbx.model import Detector

        geometry = get_root_geometry(geo=self.geo)
        shape = shape3d(self.geo)
        seg = self.geo.get_seg_geo().algo
        arows, acols = seg.asic_rows_cols()
        pix = seg.pixel_scale_size() * 1e-3
        x, y, z = geometry.get_pixel_coords(do_tilt=True)
        x.shape = y.shape = z.shape = shape
        detector = Detector()
        for n in range(shape[0]):
            for a, (r0, c0) in enumerate(seg.asic0indices()):
                nfs, nss = asic_axes(x[n], y[n], z[n], r0, c0, arows, acols)
                center = np.array((x[n,r0,c0], y[n,r0,c0], z[n,r0,c0])) * 1e-3
                origin = center - 0.5 * pix * nfs - 0.5 * pix * nss
                panel = detector.add_panel()
                panel.set_local_frame(tuple(nfs), tuple(nss), tuple(origin))
                panel.set_name("ARRAY_D0M%dA%d" % (n, a))
                panel.set_type("SENSOR_PAD")
                panel.set_pixel_size((pix, pix))
                panel.set_image_size((acols, arows))
                panel.set_trusted_range((-10.0, 2e6))
        return detector

class PsanaToCCTBX:
    """
    Class to convert psana .data geometry files to cctbx-convention dxtbx Detector models

    Emits the geometry the way cctbx.xfel does, so the result has the structure of a refined
    .expt file: one panel per psana segment, arranged in a three level hierarchy
    (detector -> module groups -> panels) whose local frames compose to the absolute positions.
    The geometry is kept in the psana laboratory frame, so the sample-detector distance from the
    input file survives and the model has a usable beam center.

    Panels follow the cctbx pixel layout, which pads each segment over the double-width pixels
    at its inter-ASIC boundaries (a 512x1024 Jungfrau segment becomes a 514x1030 panel). Raw
    psana frames must therefore be passed through ``serialtbx.detector.jungfrau.correct_panel``
    before they can be indexed with this model. Segments whose gaps are not an integral number
    of pixel pitches cannot be padded that way and are emitted one panel per ASIC instead; see
    ``panel_layout``.

    ``convert_to_file`` writes an .expt file, since that is what a cctbx consumer wants;
    ``convert`` returns the Detector object itself.

    Parameters
    ----------
    in_file : str
        Path to the input psana .data geometry file
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
    """
    def __init__(self, in_file, detname):
        self.geo = GeometryAccess(path=in_file, pbits=0, use_wide_pix_center=False)
        check_detname(detname=detname, shape=shape3d(self.geo))

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
            dxtbx Detector with a cctbx-style hierarchy
        """
        detector = cls(in_file, detname).build_detector()
        print(f"Converted psana .data geometry file {in_file} to cctbx-convention dxtbx Detector.", flush=True)
        return detector

    @classmethod
    def convert_to_file(cls, in_file, detname, out_file):
        """
        Class method to convert psana .data geometry file to a cctbx-convention .expt file

        Parameters
        ----------
        in_file : str
            Path to the input psana .data geometry file
        detname : str
            Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
        out_file : str
            Path to the output dxtbx .expt/.json geometry file

        Returns
        -------
        detector : dxtbx.model.Detector
            The Detector that was written
        """
        detector = cls.convert(in_file=in_file, detname=detname)
        write_experiment(detector=detector, out_file=out_file)
        return detector

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
        set_local_frame(node=group, basis=root_basis)
        group.set_name("D0")

        for n, module in enumerate(root.get_list_of_children()):
            subgroup = group.add_group()
            set_local_frame(node=subgroup, basis=basis_from_geo(module))
            subgroup.set_name("D0M%d" % n)

            seg_geo = sgs.Create(segname=module.oname)
            xx, yy = seg_geo.get_seg_xy_maps_um()
            xx = xx * 1e-3
            yy = yy * 1e-3
            arows, acols, indices = panel_layout(seg_geo)
            for a, (r0, c0) in enumerate(indices):
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
        self.params = load_poni(in_file)
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
        Write a psana .data file from the corrected PyFAI coordinates

        The segment records always hold image frame coordinates; ``image_frame`` only selects
        whether the IP-detector node carries the 90 degree rotation that takes them back to the
        psana laboratory frame.

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
            nfs, nss = asic_axes(x[p], y[p], z[p], 0, 0,
                                 ss_size * asics_shape[0], fs_size * asics_shape[1])
            # psana records the segment center relative to the IP-detector node, which carries
            # the mean distance
            vcent = (np.mean(x[p]), np.mean(y[p]), np.mean(z[p]) - distance)
            angle_deg_z = degrees(atan2(nfs[1], nfs[0]))
            angle_z, tilt_z = angle_and_tilt(angle_deg_z)
            tilt_x, tilt_y = tilt_xy(nfs, nss)
            angle_x = child.get_list_of_children()[p].rot_x
            angle_y = child.get_list_of_children()[p].rot_y
            recs += '\n%12s  0 %12s %2d' %(childname, self.detector.segname, p)\
                +'  %8d %8d %8d %7.0f %6.0f %6.0f   %8.5f  %8.5f  %8.5f'%\
                (vcent[0], vcent[1], vcent[2], angle_z, angle_y, angle_x, tilt_z, tilt_y, tilt_x)
        # In the laboratory frame the IP-detector node carries the 90 degree rotation
        rot_z = 0 if image_frame else 90
        recs += '\n%12s  0 %12s  0' %(topname, childname)\
            +'         0        0 %8d %7d      0      0    0.00000   0.00000   0.00000' % (distance, rot_z)
        Path(out_file).write_text(recs)

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
        self.params = load_poni(in_file)
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
        image_frame : bool
            If True, use image frame coordinates; otherwise, use psana laboratory frame coordinates
        """
        converter = cls(in_file, detector)
        converter.convert_pyfai_to_geom(out_file=out_file, image_frame=image_frame)
        print(f"Converted PyFAI .poni file {in_file} to CrystFEL .geom geometry file {out_file}.")

    def convert_pyfai_to_geom(self, out_file, image_frame=True):
        """
        Write a CrystFEL .geom file from the corrected PyFAI coordinates

        Parameters
        ----------
        out_file : str
            Path to the output .geom file
        image_frame : bool
            If True, use image frame coordinates; otherwise, use psana laboratory frame coordinates
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
        _, nasicsf = seg.number_of_asics_in_rows_cols()
        txt = header_crystfel(self.params[0])
        for n in range(nsegs):
            txt += '\n'
            for a,(r0,c0) in enumerate(seg.asic0indices()):
                nfs, nss = asic_axes(x[n], y[n], z[n], r0, c0, arows, acols)
                txt += panel_crystfel(n, a, nfs, nss, x[n,r0,c0], y[n,r0,c0], z[n,r0,c0],
                                      pix_size, arows, acols, srows, nasicsf)
                txt += tailer_crystfel(self.detector.detname)
        Path(out_file).write_text(txt)

class CrystFELToPsana:
    """
    Class to convert CrystFEL .geom geometry files to psana .data geometry files

    A .geom file carries no object tree, so the node and segment names of the psana output are
    taken from the default geometry template of the detector.

    Parameters
    ----------
    in_file : str
        Path to the CrystFEL .geom file
    detname : str
        Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
    """
    def __init__(self, in_file, detname):
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
        """
        Parse one 'key = value' line of a CrystFEL .geom file into self.dict_of_pars

        A bare key becomes a top level entry; a 'panel/key' pair becomes an entry in the
        per-panel sub-dictionary. Values are typed by key name, and lines that do not look like
        an assignment are ignored.

        Parameters
        ----------
        line : str
            A single non-comment, non-empty line of the .geom file
        """
        fields = line.split()
        nfields = len(fields)
        if nfields < 2 or fields[1] != '=':
            return
        keys = fields[0].split('/')
        nkeys = len(keys)
        if nkeys==1:
            if nfields>3:
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

    def load_geom(self, in_file):
        """
        Read a CrystFEL .geom file into the parameter dictionary self.dict_of_pars

        Parameters
        ----------
        in_file : str
            Path to the CrystFEL .geom file
        """
        self.dict_of_pars = {}
        for line in Path(in_file).read_text().splitlines():
            if not line.strip() or line.lstrip().startswith(';'):
                continue
            self._parse_line_as_parameter(line)

    def geom_to_data(self, panelasics, detname, out_file):
        """
        Write a psana .data file from the parsed CrystFEL parameters

        One record is written per panel listed in ``panelasics``, followed by the IP-detector
        record carrying the mean sample-detector distance.

        Parameters
        ----------
        panelasics : str
            Comma-separated CrystFEL panel names, one per psana segment, in segment order
        detname : str
            Detector name in pyFAI (e.g., 'Jungfrau4M', 'ePix10k2M', etc.), used for the header
        out_file : str
            Path to the output psana .data geometry file
        """
        geo = self.geo
        top = geo.get_top_geo()
        child = top.get_list_of_children()[0]
        topname = top.oname
        childname = child.oname
        seg = geo.get_seg_geo()
        sg = seg.algo
        segname = seg.oname
        X, Y, _ = sg.pixel_coord_array()
        PIX_SIZE_UM = sg.get_pix_size_um()
        M_TO_UM = 1e6
        # Offset from the first pixel of a segment to its center, from the segment template
        xc0, yc0 = X[0,0], Y[0,0]
        panels = panelasics.split(',')
        # CrystFEL stores a per-panel coffset; psana puts their mean on the IP-detector node
        segz = np.array([self.dict_of_pars[panel].get('coffset', 0) for panel in panels])
        distance = self.dict_of_pars.get('coffset', 0) + round(segz.mean()*1e6)*1e-6
        recs = header_psana(detname=detname)
        for p, panel in enumerate(panels):
            dicasic = self.dict_of_pars[panel]
            nfs = np.array(dicasic['fs'], dtype=np.float64)
            nss = np.array(dicasic['ss'], dtype=np.float64)
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
        Path(out_file).write_text(recs)

    def convert_geom_to_data(self, detname, out_file):
        """
        Write a psana .data file from the loaded .geom file

        Parameters
        ----------
        detname : str
            Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)
        out_file : str
            Path to the output psana .data geometry file
        """
        panelasics = panel_names(detname)
        if panelasics is None:
            raise ValueError(f"Detector name {detname} not implemented.")
        self.geom_to_data(panelasics, get_pyFAI_detname(detname), out_file)

class CrystFELToPyFAI:
    """
    Class to convert CrystFEL .geom geometry files to PyFAI Detector model by using intermediate psana .data files
    """

    @staticmethod
    def convert(in_file, detname):
        """
        Convert a CrystFEL .geom geometry file to a PyFAI Detector instance

        Parameters
        ----------
        in_file : str
            Path to the CrystFEL .geom file
        detname : str
            Detector name in psana (e.g., 'jungfrau4M', 'epix10k2M', etc.)

        Returns
        -------
        detector : pyFAI.Detector
            PyFAI Detector instance with geometry defined by the input .geom file
        """
        fd, temp_file = tempfile.mkstemp(suffix=".data")
        os.close(fd)
        try:
            CrystFELToPsana.convert(in_file=in_file, detname=detname, out_file=temp_file)
            detector = PsanaToPyFAI.convert(in_file=temp_file, detname=detname)
        finally:
            os.unlink(temp_file)
        print(f"Converted CrystFEL .geom file {in_file} to PyFAI Detector instance.", flush=True)
        return detector

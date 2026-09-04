"""
Geometry correction functions for LCLS detectors based on pyFAI calibration parameters.
"""

import numpy as np
from math import degrees, atan2

def rotation_matrix(params):
    """
    Compute and return the detector tilts as a single rotation matrix

    Parameters
    ----------
    params : list
        Detector parameters found by PyFAI calibration
    """
    cos_rot1 = np.cos(params[3])
    cos_rot2 = np.cos(params[4])
    cos_rot3 = np.cos(params[5])
    sin_rot1 = np.sin(params[3])
    sin_rot2 = np.sin(params[4])
    sin_rot3 = np.sin(params[5])
    # Rotation about axis 1: Note this rotation is left-handed
    rot1 = np.array([[1.0, 0.0, 0.0],
                        [0.0, cos_rot1, sin_rot1],
                        [0.0, -sin_rot1, cos_rot1]])
    # Rotation about axis 2. Note this rotation is left-handed
    rot2 = np.array([[cos_rot2, 0.0, -sin_rot2],
                        [0.0, 1.0, 0.0],
                        [sin_rot2, 0.0, cos_rot2]])
    # Rotation about axis 3: Note this rotation is right-handed
    rot3 = np.array([[cos_rot3, -sin_rot3, 0.0],
                        [sin_rot3, cos_rot3, 0.0],
                        [0.0, 0.0, 1.0]])
    rotation_matrix = np.dot(np.dot(rot3, rot2), rot1)  # 3x3 matrix
    return rotation_matrix

def get_beam_center(params):
    """
    From pyFAI calibration parameters, return the beam center in meters

    Parameters
    ----------
    params : list
        Detector parameters found by PyFAI calibration
    """
    dist = params[0]
    poni1 = params[1]
    poni2 = params[2]
    rot1 = params[3]
    rot2 = params[4]
    distance = dist / (np.cos(rot1) * np.cos(rot2))
    cx = poni1 + dist * np.tan(rot2) / np.cos(rot1)
    cy = poni2 - dist * np.tan(rot1)
    return distance, cx, cy

def angle_and_tilt(a):
    """
    From a given angle, return the angle and tilt in degrees

    Parameters
    ----------
    a : float
        Angle in degrees
    """
    desangles = np.array((-180,-90, 0, 90, 180))
    difangles = a-desangles
    absdifang = np.abs(difangles)
    imin = np.argmin(absdifang)
    angle, tilt = desangles[imin], difangles[imin]
    return (angle if angle>=0 else angle+360), tilt

def unit_vector_pitch_angle_max_ind(u):
    """
    From a given unit vector, return the pitch angle and index of maximum value
    in the vector

    Parameters
    ----------
    u : numpy.ndarray
        Input unit vector
    """
    absu = np.abs(u)
    imax = np.where(absu == np.amax(absu))[0]
    pitch = degrees(atan2(u[2],u[imax]))
    pitch = (pitch+180) if pitch<-90 else (pitch-180) if pitch>90 else pitch
    return pitch, imax

def tilt_xy(uf, us):
    """
    From two unit vectors, return the tilt in x and y directions
    based on the maximum index

    Parameters
    ----------
    uf : numpy.ndarray
        First unit vector
    us : numpy.ndarray
        Second unit vector
    """
    tilt_f, imaxf = unit_vector_pitch_angle_max_ind(uf)
    tilt_s, imaxs = unit_vector_pitch_angle_max_ind(us)
    tilt_x, tilt_y = (tilt_s, tilt_f) if imaxf == 0 else (tilt_f, tilt_s)
    return tilt_x, -tilt_y

def rotate_z(angle_z, angle_x, angle_y, tilt_x, tilt_y):
    """
    For a given angle around Z-axis, switch angle_x and angle_y appropriately.

    Parameters
    ----------
    angle_z : float
        Angle around Z-axis in degrees in [0, 90, 180, 270]
    angle_x : float
        Angle around X-axis in degrees in [0, 90, 180, 270]
    angle_y : float
        Angle around Y-axis in degrees in [0, 90, 180, 270]
    tilt_x : float
        Tilt angle around X-axis in degrees
    tilt_y : float
        Tilt angle around Y-axis in degrees
    """
    if angle_z == 0.0:
        result_x, result_y = angle_x, angle_y
    elif angle_z == 90.0:
        result_x, result_y = angle_y, -angle_x
        tilt_x, tilt_y = tilt_y, -tilt_x
        result_y = result_y if result_y >= 0 else result_y + 360
    elif angle_z == 180.0:
        result_x, result_y = -angle_x, -angle_y
        tilt_x, tilt_y = -tilt_x, -tilt_y
        result_x = result_x if result_x >= 0 else result_x + 360
        result_y = result_y if result_y >= 0 else result_y + 360
    elif angle_z == 270.0:
        result_x, result_y = -angle_y, angle_x
        tilt_x, tilt_y = -tilt_y, tilt_x
        result_x = result_x if result_x >= 0 else result_x + 360
    else:
        raise ValueError("angle_z must be one of [0, 90, 180, 270]")
    result_x = 0.0 if result_x == 0.0 else result_x
    result_y = 0.0 if result_y == 0.0 else result_y
    return result_x, result_y, tilt_x, tilt_y

def correct_geometry(detector, params):
    """
    Correct the geometry based on the PONI parameters

    Parameters
    ----------
    detector : pyFAI.detectors.Detector
        The pyFAI Detector object to be corrected
    params : list
        Detector parameters found by PyFAI calibration
    """
    p1, p2, p3 = detector.calc_cartesian_positions()
    dist = params[0]
    poni1 = params[1]
    poni2 = params[2]
    p1 = (p1 - poni1).ravel()
    p2 = (p2 - poni2).ravel()
    if p3 is None:
        p3 = np.zeros_like(p1) + dist
    else:
        p3 = (p3+dist).ravel()
    coord_det = np.stack((p1, p2, p3), axis=0)
    coord_sample = np.dot(rotation_matrix(params), coord_det)
    x, y, z = coord_sample
    return x, y, z

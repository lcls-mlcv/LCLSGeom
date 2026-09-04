"""
Formatting utilities for LCLS geometry files.
"""

import getpass
from time import strftime, localtime
from LCLSGeom.utils.calibmap import detname_to_groups, detname_to_quadrants

def str_to_int_or_float(s):
    """
    Checks if s is a float or int and converts to corresponding type.

    Parameters
    ----------
    s : str
        Input string to check
    """
    v = float(s)
    if v%1 == 0: v=int(v)
    return v

def sfields_to_xyz_vector(flds):
    """ 
    Converts ['+0.002583x', '-0.999997y', '+0.000000z'] to (0.002583, -0.999997, 0.000000)

    Parameters
    ----------
    flds : list
        List of strings representing x, y, and z coordinates
    """
    v = (float(flds[0].strip('x')), float(flds[1].strip('y')))
    z = float(flds[2].strip('z')) if len(flds)==3 else 0
    v += (z,)
    return v

def str_is_segment_and_asic(s):
    """ 
    Checks if s looks like str 'q0a2' or 'p12a7' and
    returns 'p0.2' or 'p12.7' or False

    Parameters
    ----------
    s : str
        Input string to check
    """
    if not isinstance(s, str)\
    or len(s)<2: return False
    flds = s[1:].split('a')
    return False if len(flds) !=2 else\
        'p%sa%s' % (flds[0], flds[1]) if all([f.isdigit() for f in flds]) else\
        False

def get_time_stamp(fmt='%Y-%m-%dT%H:%M:%S', time=None):
    """
    Returns the current timestamp formatted according to the specified format.

    Parameters
    ----------
    fmt : str, optional
        The format string for the timestamp (default is '%Y-%m-%dT%H:%M:%S')
    time : float, optional
        The time in seconds since the epoch (default is None, which uses the current time)
    """
    return strftime(fmt, localtime(time))

def get_login():
    """
    Returns the login name of the current user.

    Parameters
    ----------
    None
    """
    return getpass.getuser()

def header_psana(detname):
    """
    Returns the header for the psana geometry file based on the detector name.
    The header includes information such as the title, date, author, experiment,
    detector, calibration type, and comments.

    Parameters
    ----------
    detname : str
        The name of the detector (e.g., 'Rayonix', 'ePix10k2M', etc.)
    """
    if detname.lower() == 'rayonix':
        txt=\
        '# TITLE      Geometry parameters of Rayonix'\
        +'\n# DATE_TIME  %s' % get_time_stamp(fmt='%Y-%m-%dT%H:%M:%S %Z')\
        +'\n# METROLOGY  no metrology available'\
        +'\n# AUTHOR     %s' % get_login()\
        +'\n# EXPERIMENT N/A'\
        +'\n# DETECTOR   Rayonix'\
        +'\n# CALIB_TYPE geometry'\
        +'\n# COMMENT:01 Automatically created from BayFAI for the Rayonix detector'\
        +'\n'
    elif detname.lower() == 'epix10k2m':
        txt=\
        '# TITLE      Geometry parameters of ePix10k2M'\
        +'\n# DATE_TIME  %s' % get_time_stamp(fmt='%Y-%m-%dT%H:%M:%S %Z')\
        +'\n# METROLOGY  no metrology available'\
        +'\n# AUTHOR     %s' % get_login()\
        +'\n# EXPERIMENT N/A'\
        +'\n# DETECTOR   Epix10ka2M'\
        +'\n# CALIB_TYPE geometry'\
        +'\n# COMMENT:01 Automatically created from BayFAI for the 16-segment ePix10k2M detector'\
        +'\n'
    elif detname.lower() == 'epix10kaquad':
        txt=\
        '# TITLE      Geometry parameters of ePix10kaQuad'\
        +'\n# DATE_TIME  %s' % get_time_stamp(fmt='%Y-%m-%dT%H:%M:%S %Z')\
        +'\n# METROLOGY  no metrology available'\
        +'\n# AUTHOR     %s' % get_login()\
        +'\n# EXPERIMENT N/A'\
        +'\n# DETECTOR   %s' % detname\
        +'\n# CALIB_TYPE geometry'\
        +'\n# COMMENT:01 Automatically created from BayFAI for the ePix10kaQuad detector'\
        +'\n'
    elif detname.lower() == 'jungfrau05m':
        txt=\
        '# TITLE      Geometry parameters of Jungfrau1M'\
        +'\n# DATE_TIME  %s' % get_time_stamp(fmt='%Y-%m-%dT%H:%M:%S %Z')\
        +'\n# METROLOGY  no metrology available'\
        +'\n# AUTHOR     %s' % get_login()\
            +'\n# EXPERIMENT N/A'\
        +'\n# DETECTOR   DetLab.0:Jungfrau.0 or jungfrau05M'\
        +'\n# CALIB_TYPE geometry'\
        +'\n# COMMENT:01 Automatically created from BayFAI for the 1-segment Jungfrau1M detector'\
        +'\n'
    elif detname.lower() == 'jungfrau1m':
        txt=\
        '# TITLE      Geometry parameters of Jungfrau1M'\
        +'\n# DATE_TIME  %s' % get_time_stamp(fmt='%Y-%m-%dT%H:%M:%S %Z')\
        +'\n# METROLOGY  no metrology available'\
        +'\n# AUTHOR     %s' % get_login()\
        +'\n# EXPERIMENT N/A'\
        +'\n# DETECTOR   DetLab.0:Jungfrau.1 or jungfrau1M'\
        +'\n# CALIB_TYPE geometry'\
        +'\n# COMMENT:01 Automatically created from BayFAI for the 2-segment Jungfrau1M detector'\
        +'\n'
    elif detname.lower() == 'jungfrau4m':
        txt=\
        '# TITLE      Geometry parameters of Jungfrau4M'\
        +'\n# DATE_TIME  %s' % get_time_stamp(fmt='%Y-%m-%dT%H:%M:%S %Z')\
        +'\n# METROLOGY  no metrology available'\
        +'\n# AUTHOR     %s' % get_login()\
        +'\n# EXPERIMENT N/A'\
        +'\n# DETECTOR   DetLab.0:Jungfrau.2 or jungfrau4M'\
        +'\n# CALIB_TYPE geometry'\
        +'\n# COMMENT:01 Automatically created from BayFAI for the 8-segment Jungfrau4M detector'\
        +'\n'
    elif detname.lower() == 'jungfrau16m':
        txt=\
        '# TITLE      Geometry parameters of Jungfrau16M'\
        +'\n# DATE_TIME  %s' % get_time_stamp(fmt='%Y-%m-%dT%H:%M:%S %Z')\
        +'\n# METROLOGY  no metrology available'\
        +'\n# AUTHOR     %s' % get_login()\
        +'\n# EXPERIMENT N/A'\
        +'\n# DETECTOR   DetLab.0:Jungfrau.3 or jungfrau16M'\
        +'\n# CALIB_TYPE geometry'\
        +'\n# COMMENT:01 Automatically created from BayFAI for the 32-segment Jungfrau16M detector'\
        +'\n'
    else:
        txt=''
    txt +=\
        '\n# PARAM:01 PARENT     - name and version of the parent object'\
        '\n# PARAM:02 PARENT_IND - index of the parent object'\
        '\n# PARAM:03 OBJECT     - name and version of the object'\
        '\n# PARAM:04 OBJECT_IND - index of the new object'\
        '\n# PARAM:05 X0         - x-coordinate [um] of the object origin in the parent frame'\
        '\n# PARAM:06 Y0         - y-coordinate [um] of the object origin in the parent frame'\
        '\n# PARAM:07 Z0         - z-coordinate [um] of the object origin in the parent frame'\
        '\n# PARAM:08 ROT_Z      - object design rotation angle [deg] around Z axis of the parent frame'\
        '\n# PARAM:09 ROT_Y      - object design rotation angle [deg] around Y axis of the parent frame'\
        '\n# PARAM:10 ROT_X      - object design rotation angle [deg] around X axis of the parent frame'\
        '\n# PARAM:11 TILT_Z     - object tilt angle [deg] around Z axis of the parent frame'\
        '\n# PARAM:12 TILT_Y     - object tilt angle [deg] around Y axis of the parent frame'\
        '\n# PARAM:13 TILT_X     - object tilt angle [deg] around X axis of the parent frame'\
        '\n'\
        '\n# HDR PARENT IND     OBJECT IND    X0[um]   Y0[um]   Z0[um]   ROT-Z  ROT-Y  ROT-X     TILT-Z    TILT-Y    TILT-X'
    return txt

def header_crystfel(distance):
    """
    Returns the header for the CrystFEL geometry file.
    The header includes the clen and photon energy PVs.

    Parameters    
    ----------
    distance : float
        The distance from the source to the detector [m]
    """
    txt = (
        f'\n; Geometry Constants generated by LCLSGeom'
        f'\n'
        f'\nclen = {distance:f}'
        f'\nphoton_energy = /LCLS/photon_energy_eV'
        f'\nadu_per_eV = 0.1'
        f'\n'
        f'\ndata = /entry_1/data_1/data'
        f'\n'
        f'\ndim0 = %'
        f'\ndim1 = ss'
        f'\ndim2 = fs'
        f'\n'
        f'\n; mask = /entry_1/data_1/mask'
        f'\n; mask_good = 0x0000'
        f'\n; mask_bad = 0xffff'
    )
    return txt

def tailer_crystfel(detname):
    """
    Returns the tailer for the CrystFEL geometry file.
    The tailer includes group definitions for millipede.

    Parameters
    ----------
    detname : str
        The name of the detector in pyFAI (e.g., 'Jungfrau16M', 'ePix10k2M', etc.)
    """
    groups = detname_to_groups.get(detname.lower())
    quadrants = detname_to_quadrants.get(detname.lower())

    lines = []

    if groups:
        for group_name, members in groups.items():
            lines.append(f'{group_name} = {members}')

    if quadrants:
        for group_name, members in quadrants.items():
            lines.append(f'{group_name} = {members}')

    return '\n' + '\n'.join(lines)
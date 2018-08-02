#!/usr/bin/env python3
# author: Joseph Lucero
# created on: 28 June 2018 17:49:24
# purpose: functions that acquire information 

from __future__ import division

from sys import exit
from os import getcwd, environ
from datetime import date
from os.path import isdir, exists,expanduser

from numpy import loadtxt, genfromtxt, array
from pandas import read_csv

from pydicom import dcmread

def get_target_dir():

    """\
    Description:

    Inputs:

    Outputs:
    """

    while True:
        target_dir = expanduser(input(
        "Please input path to directory in which all relevant files are " + \
        "stored :>> """
        ))
        if not isdir(target_dir):
            print("That is not a valid path to a directory.")
            answer = input(
                "Should we set the current working directory as path? [y/n] :>>"
                )
            if answer.lower() in ('y', 'yes'):
                target_dir = getcwd()
                break
            elif answer.lower() in ('n', 'no'):
                continue
            else:
                print("Input not understood.")
                continue
        else:
            break

    return target_dir

def get_ref_info_from_ref_slice(ref_file_1, ref_file_2, num_CT_files):

    """\
    Description:
    Acquire the necessary reference info that is to be used to create the 
    egsphant from the CT image

    Inputs:
    :param ref_file_1: First reference file path (default lowest inferior CT file) 
    :type ref_file_1: string
    :param ref_file_2: Second reference file path (default second lowest inferior CT file)
    :type ref_file_2: string
    :param num_CT_files: Total number of CT files
    :type num_CT_files: int

    Outputs:
    :param reference_dict: dictionary of reference information
        "SLICE THICKNESS" : Distance between each CT image
        "START OF GRID" : The lowest leftmost point in the most inferior CT image
        "SIZE OF GRID" : Number of voxels in x, y, and z
        "VOXEL DIMS" : Dimensions of voxels
        "INTERCEPT" : Intercept that is applied to the CT images reshift to appropriate values
        "SCALING" : ??? (Not sure what this does)
        "BOUNDS" : Coordinates of voxel boundaries
        "VOXEL CENTERS" : Location of the voxel ceneters 
    :type reference_dict: dict\
    """
    ref_1 = dcmread(ref_file_1, force=True)
    ref_2 = dcmread(ref_file_2, force=True)

    assert ref_1.RescaleIntercept == ref_2.RescaleIntercept, \
        "RED ALERT: Intercepts are different for different reference files!!"

    SLICE_THICKNESS = float(ref_2.SliceLocation - ref_1.SliceLocation)
    START_OF_GRID = tuple(
        map(
            float,
            (
                ref_1.ImagePositionPatient[0],
                ref_1.ImagePositionPatient[1],
                ref_1.SliceLocation
            )
        )
    )
    SIZE_OF_GRID = tuple(
        map(
            int,
            (
                ref_1.Columns, ref_1.Rows, num_CT_files
            )
        )
    )
    VOXEL_DIMS = tuple(
        map(
            float,
            (
                ref_1.PixelSpacing[0],
                ref_1.PixelSpacing[1],
                SLICE_THICKNESS
            )
        )
    )
    INTERCEPT = float(ref_1.RescaleIntercept)
    SCALING = float(ref_1.RescaleSlope)

    print("Size of grid:", SIZE_OF_GRID)
    print("Voxel dimensions:", VOXEL_DIMS)

    # Get voxel bounds
    XBOUNDS, YBOUNDS, ZBOUNDS = [], [], []  # initialize the lists
    # set the boundary values
    XBOUNDS.append((START_OF_GRID[0] - (0.5 * VOXEL_DIMS[0])) / 10)
    YBOUNDS.append((START_OF_GRID[1] - (0.5 * VOXEL_DIMS[1])) / 10)
    ZBOUNDS.append((START_OF_GRID[2] - (0.5 * VOXEL_DIMS[2])) / 10)

    for i in range(SIZE_OF_GRID[0]):
        XBOUNDS.append(XBOUNDS[i] + (VOXEL_DIMS[0] / 10))
    for j in range(SIZE_OF_GRID[1]):
        YBOUNDS.append(YBOUNDS[j] + (VOXEL_DIMS[1] / 10))
    for k in range(SIZE_OF_GRID[2]):
        ZBOUNDS.append(ZBOUNDS[k] + (VOXEL_DIMS[2] / 10))

    VOXEL_CENTERS = []
    for y_index in range(len(YBOUNDS) - 1):
        for x_index in range(len(XBOUNDS) - 1):
            VOXEL_CENTERS.append(
                (
                    XBOUNDS[x_index] + ((0.5 * VOXEL_DIMS[0]) / 10),
                    YBOUNDS[y_index] + ((0.5 * VOXEL_DIMS[1]) / 10)
                )
            )

    return {
        "SLICE THICKNESS": SLICE_THICKNESS, "START OF GRID": START_OF_GRID,
        "SIZE OF GRID": SIZE_OF_GRID, "VOXEL DIMS": VOXEL_DIMS,
        "INTERCEPT": INTERCEPT, "SCALING": SCALING,
        "BOUNDS": (XBOUNDS, YBOUNDS, ZBOUNDS),
        "VOXEL CENTERS": VOXEL_CENTERS
    }

def get_CT_calibration(path_to_calibration=None):
    
    """\
    Description:
    
    Inputs:

    Outputs:\
    """

    cwd = getcwd()  # This is not wise. Need to figure out what the best way to call this is
    
    if not path_to_calibration:
        while True:
            path_to_calibration = input(
                "Please input the FULL file path to the CT calibration data:>> "
            )
            if not exists(path_to_calibration):
                print(
                    "This is not a valid file path to an existing " 
                    + "calibration file."
                    )
                answer = input("Use default (Peppa et al.)? [y/n]:>> ")
                if answer.lower() in ('y','yes'):
                    path_to_calibration = cwd + '/lib/calibration/PeppaCalibrationCurve.dat'
                elif answer.lower() in ('n', 'no'):
                    pass
                else: 
                    print("Input not understood. Try again.")
            else:
                break
    else:
        if path_to_calibration.lower() == 'default':
            print("Calling default calibration curve.")
            path_to_calibration = cwd + '/lib/calibration/PeppaCalibrationCurve2.dat'

    HU, mass_density = loadtxt(path_to_calibration,unpack=True)

    return HU, mass_density

def get_media(path_to_media=None):

    """\
    Description:

    Inputs:

    Outputs:\
    """

    cwd = getcwd()  # This is not wise. Need to figure out what the best way to call this is

    if not path_to_media:
        while True:
            path_to_media = input(
                "Please input the FULL file path to the media definitions file curve:>> "
            )
            if not exists(path_to_media):
                print(
                    "This is not a valid file path to an existing "
                    + "calibration file."
                    )
                answer = input("Use default (Peppa et al.)? [y/n]:>> ")
                if answer.lower() in ('y', 'yes'):
                    path_to_media = cwd + '/lib/calibration/PeppaMedia.dat'
                elif answer.lower() in ('n', 'no'):
                    pass
                else:
                    print("Input not understood. Try again.")
            else:
                break
    elif path_to_media.lower() == 'default':
        print("Calling default media file.")
        path_to_media = cwd + '/lib/media/PeppaMedia.dat'

    media_info = read_csv(path_to_media, delim_whitespace=True, header=None)

    media_name = list(media_info[0])
    media_HU_min = array(media_info[1])
    media_HU_max = array(media_info[2])
    media_HU = array(media_info[3])
    media_density = array(media_info[4])
    PEGS4_name = list(media_info[5])

    media_HU_ranges = [
        (HU_min, HU_max) for HU_min, HU_max in zip(media_HU_min, media_HU_max)
    ]

    return media_name, media_HU_ranges, media_HU, media_density, PEGS4_name

def get_egs_brachy_settings(file_source,**kwargs):

    """\
    Description: Pull the settings needed for the egs_brachy simulation from 
    file input by user

    Inputs:
    :param :
    :type :
    
    Outputs:
    :param :
    :type :
    \
    """

    # Call environment variable EGS_BRACHY pointing to egs_brachy directory
    egs_brachy = environ['EGS_BRACHY']  

    format_dict = {
        'date': date.today(),
        'info_dir': getcwd(),
        'file_source': file_source,
    }

    input_settings = {
        'run_mode': 'normal',

        'AE': 1.512,
        'UE': 2.012,
        'AP': 0.001,
        'UP': 1.500,

        'material_dat_file_path': egs_brachy + '/lib/media/material.dat',
        'egsphant_file_path': 
        egs_brachy + '/lib/geometry/phantoms/egsphant/PeppaBreastHDR192Ir_MBDCA-WG.egsphant.gz',
        'density_file_path': egs_brachy + '/lib/media/material.dat',
        'seed_file_path': 
        egs_brachy + '/lib/geometry/sources/Ir192_HDR/MBDCA-WG/MBDCA-WG.geom',
        'transform_file_path':
        egs_brachy + '/lib/geometry/transformations/PeppaBreastHDR192IR_MBDCA-WGG_srcPosnRotn', 

        'volcor_density': 1E8,
        'seed_boundary_file_path':
        egs_brachy + '/lib/geometry/sources/Ir192_HDR/MBDCA-WG/boundary.shape',
        'volcor_type': 'none',

        'seed_shape_file_path': 
        egs_brachy + '/lib/geometry/sources/Ir192_HDR/MBDCA-WG/MBDCA-WG.shape',
        'spectrum_file_path': 
        egs_brachy + '/lib/spectra/Ir192_NNDC_2.6_line.spectrum',

        'score_tracklength': 'yes',
        'score_deposition': 'no',
        'score_scatter': 'no',
        'muen_file_path': egs_brachy + '/lib/media/material.dat',
        'output_egsphant': 'no',
        'output_voxel_info': 'no',
        'output_volcor_file': 'no',
        'record_init_pos': 0.0,
        'dose scaling_factor': 1.0,

        'transport_param_file_path': 
        egs_brachy + '/lib/transport/high_energy_default',
        }

    return format_dict, input_settings




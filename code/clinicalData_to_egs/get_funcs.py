#!/usr/bin/env
# author: Joseph Lucero
# created on: 28 June 2018 17:49:24
# purpose: functions that acquire information 

from __future__ import division

from os import getcwd
from os.path import isdir, exists

from numpy import loadtxt, genfromtxt

from dicom import read_file


def get_defaults():
    """
    Description:
    Define some default settings that are needed for the code.
    Return settings in a dict that is accessible to the user

    Inputs:
    None

    Outputs:
    """

    NAME_STRING = "AF_FINAL"
    NUM_OF_HISTORIES = "1E9"
    AG_CUTOFF_DEN = 0.985
    CALC_CUTOFF_DEN = 1.16
    BONE_CUTOFF_DEN = 1.16
    LOW_ARTIFACT_CUTOFF_DEN = 1.066
    CARTILAGE_CUTOFF_DEN = 1.066
    YELLOW_MARROW_CUTOFF_DEN = 1.005
    RED_MARROW_CUTOFF_DEN = 1.055
    INDIVIDUALIZATION = False

    return {
        "NAME_STRING": NAME_STRING, "NUM_OF_HISTORIES": NUM_OF_HISTORIES,
        "AG_CUTOFF_DEN": AG_CUTOFF_DEN, "CALC_CUTOFF_DEN": CALC_CUTOFF_DEN,
        "BONE_CUTOFF_DEN": BONE_CUTOFF_DEN,
        "LOW_ARTIFACT_CUTOFF_DEN": LOW_ARTIFACT_CUTOFF_DEN,
        "CARTILAGE_CUTOFF_DEN": CARTILAGE_CUTOFF_DEN,
        "YELLOW_MARROW_CUTOFF_DEN": YELLOW_MARROW_CUTOFF_DEN,
        "RED_MARROW_CUTOFF_DEN": RED_MARROW_CUTOFF_DEN,
        "INDIVIDUALIZATION": INDIVIDUALIZATION
    }


def get_target_dir():
    while True:
        target_dir = raw_input(
        "Please input path to directory in which all relevant files are " + \
        "stored :>> """
        )
        if not isdir(target_dir):
            print "That is not a valid path to a directory."
            answer = raw_input(
                "Should we set the current working directory as path? [y/n] :>>"
                )
            if answer.lower() in ('y', 'yes'):
                target_dir = getcwd()
                break
            elif answer.lower() in ('n', 'no'):
                continue
            else:
                print "Input not understood."
                continue
        else:
            break

    return target_dir

def get_ref_info_from_ref_slice(ref_file_1, ref_file_2, num_CT_files):

    """
    Description:

    Inputs:

    Outputs:
    """
    ref_1 = read_file(ref_file_1)
    ref_2 = read_file(ref_file_2)
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

    print "Voxel dimensions:", VOXEL_DIMS

    # Get voxel bounds
    XBOUNDS, YBOUNDS, ZBOUNDS = [], [], []  # initialize the lists
    # set the boundary values
    XBOUNDS.append((START_OF_GRID[0] - ((0.5 * VOXEL_DIMS[0])/10)))
    YBOUNDS.append((START_OF_GRID[1] - ((0.5 * VOXEL_DIMS[1])/10)))
    ZBOUNDS.append((START_OF_GRID[2] - ((0.5 * VOXEL_DIMS[2])/10)))

    for i, j, k in zip(
            xrange(SIZE_OF_GRID[0]),
            xrange(SIZE_OF_GRID[1]),
            xrange(SIZE_OF_GRID[2])):

        XBOUNDS.append(XBOUNDS[i] + (VOXEL_DIMS[0] / 10))
        YBOUNDS.append(YBOUNDS[j] + (VOXEL_DIMS[1] / 10))
        ZBOUNDS.append(ZBOUNDS[k] + (VOXEL_DIMS[2] / 10))

    VOXEL_CENTERS = []
    for y_index in xrange(len(YBOUNDS) - 1):
        for x_index in xrange(len(XBOUNDS) - 1):
            VOXEL_CENTERS.append(
                (
                    XBOUNDS[x_index] + ((0.5 * VOXEL_DIMS[0]) / 10),
                    YBOUNDS[y_index] + ((0.5 * VOXEL_DIMS[1]) / 10)
                )
            )

    return {
        "SLICE_THICKNESS": SLICE_THICKNESS, "START_OF_GRID": START_OF_GRID,
        "SIZE_OF_GRID": SIZE_OF_GRID, "VOXEL_DIMS": VOXEL_DIMS,
        "INTERCEPT": INTERCEPT, "SCALING": SCALING,
        "BOUNDS": (XBOUNDS, YBOUNDS, ZBOUNDS),
        "VOXEL_CENTERS": VOXEL_CENTERS
    }

def get_CT_calibration(path_to_calibration=None):
    """
    Description:
    
    Inputs:

    Outputs:
    """
    
    if not path_to_calibration:
        while True:
            path_to_calibration = raw_input(
                "Please input the FULL file path to the CT calibration curve:>> "
            )
            if not exists(path_to_calibration):
                print \
                "This is not a valid file path to an existing calibration file." + \
                " Please specify a correct path."
            else:
                break

    HU, mass_density = loadtxt(path_to_calibration,unpack=True)

    return HU, mass_density

def get_media(path_to_media=None):

    if not path_to_media:
        while True:
            path_to_media = raw_input(
                "Please input the FULL file path to the media definitions file curve:>> "
            )
            if not exists(path_to_media):
                print \
                    "This is not a valid file path to an existing calibration file." + \
                    " Please specify a correct path."
            else:
                break

    [
        media_name, media_HU, media_density, PEGS4_name, HU_min, HU_max
        ] = genfromtxt(path_to_media, unpack=True)

    return media_name, media_HU, media_density, PEGS4_name, HU_min, HU_max




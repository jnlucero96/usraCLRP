#!/usr/bin/env python3
# author: Joseph Lucero
# created on: 28 June 2018 17:43:28
# purpose: functions that process files and extract information

# import python modules
from __future__ import division

from os import walk
from glob import glob
from pickle import load as pload
from os.path import join, isfile, expanduser
from pydicom import dcmread

def process_CT(root_dir):

    """\
    Description: 
    Opening sorting, and arraying the CT files.
    
    Inputs:
    :param root_dir: directory in which all CT files are stored
    :type root_dir: string

    Outputs:
    :param ordered_CT_list: ordered list of dicom files that are ordered from \
    inferior to superior
    :type ordered_CT_list: list
    :param slice_array_lst: every entry corresponds to the intensity data \
    for the corresponding CT image
    :type slice_array_lst: list of arrays
    :param n_CT: number of CT files that are in the root directory
    :type n_CT: int\
    """

    print('\n\n' + '=' * 50)
    print("Now importing and arranging CT data")
    print('=' * 50)

    # extract the CT filenames in the target directory
    filename_CT_lst = glob(join(root_dir, '*CT*.dcm'))

    # for CT_file in filename_CT_lst:
    #     print CT_file, dcmread(CT_file).SliceLocation

    sorted_dcm = sorted([
        (
            dcmread(CT_file).SliceLocation,
            CT_file
        ) for CT_file in filename_CT_lst 
    ])  # load in the files

    # sort the dcm files from inferior (-) to superior (+) and
    # extract the ordered CT Dicom files    
    ordered_CT_lst = [
        item[1] for item in sorted_dcm
    ]

    slice_array_lst = [
        dcmread(current_CT).pixel_array for current_CT in ordered_CT_lst
    ]  # slice array list indexing goes [z][y][x]

    return ordered_CT_lst, slice_array_lst, len(ordered_CT_lst)


def process_RP(root_dir, scale_dose=False):

    """\
    Description: 
    Open the plan dicom files, process them, and extract relevant data.
    
    Inputs:
    :param root_dir: directory in which the plan dicom files are stored
    :type root_dir: str
    :param scale_dose: parameter that determines whether dose inside the plan file will be scaled or not.
    :type scale_dose: bool

    Outputs:
    :param final_seed_locations: locations of all the seeds that are within the phantom
    :type final_seed_locations: list
    :param DSF: scaling factor that turns dose output by the TPS and egs_brachy to absolute dose in Gy
    :type DSF: float64\
    """

    # extract the filenames in the RT directory
    filename_RP_lst = glob(join(root_dir, '*RP*.dcm')); answer = 0 

    num_RP_files = len(filename_RP_lst)

    # there should only be a unique RT (RP) plan file for a CT directory.
    # this covers the cases when user has more than one in the directory 
    if num_RP_files > 1:
        print("There is more than one RT plan in this directory:")
        for index, file_name in enumerate(filename_RP_lst):
            print(str(index) + ": " + file_name)
        answer = int(input(
            "Which of the {0} files do you wish to use ".format(index + 1)
            + "as the plan file? Please input in the index [int]:>> "
        ))

    current_plan_file = filename_RP_lst[answer]

    rp = dcmread(current_plan_file)
    print(
        "Now finding treatement plan info. for Patient {0}.".format(rp.PatientID)
        )

    try:
        air_kerma = float(rp.SourceSequence[0].ReferenceAirKermaRate)
        print('Air Kerma =', air_kerma)
    except (NameError, AttributeError):
        print("No Air Kerma found for Patient {0}.".format(rp.PatientID))
        scale_dose = False  # without air kerma cannot scale dose

    # define dose scaling factor (DSF)
    if scale_dose:
        # what is this hard coded number? for what source?
        # from Stephen this is for Pd-103
        DSF = 91558000000000 * air_kerma
    else:
        DSF = 1.0
        print("Dose remains unscaled.")

    try:
        # get final list of seed locations
        final_seed_locations = [
            [
                float(
                    seed_info.ChannelSequence[0].BrachyControlPointSequence[1].ControlPoint3DPosition[i]
                )
                for i in range(0, 3)
            ] for __, seed_info in enumerate(rp.ApplicationSetupSequence)
        ]
        print("Seed locations found.")
    except:
        print("No seed locations.")

    return final_seed_locations, DSF

def process_RS(root_dir, SIZE_OF_GRID):

    """\
    Description: 
    Takes the contour file and loads them into the program so that it can aid media assignment later

    Inputs:
    :param root_dir: Directory in which the contour file is stored
    :type root_dir: str
    :param SIZE_OF_GRID: Resolution of the CT images x Number of CT images 
    :type SIZE_OF_GRID: tuple of ints

    Outputs:
    :param cont_map: Contours for every valid structure that is found within the given files
    :type cont_map: list of lists\
    """

    # glob the contour files in the directory
    contour_files = glob(join(root_dir, '*RS*.dcm')); answer = 0

    num_RS_files = len(contour_files)
    
    # there should only be a unique contour (RS) file for a CT directory.
    # this covers the cases when user has more than one in the directory
    if num_RS_files > 1:
        print("There is more than one RS plan in this directory:")
        for index, file_name in enumerate(contour_files):
            print(str(index) + ": " + file_name)
        answer = int(input(
            "Which of the {0} files do you wish to use ".format(index + 1)
            + "as the plan file? Please input in the index [int]:>> "
        ))

    current_contour_file = contour_files[answer]

    rt = dcmread(current_contour_file)

    structure_names = [
        name.ROIName for name in rt.StructureSetROISequence
    ]

    if structure_names:
        print("Structures Found:")
        print(structure_names)

    ROI_DICT = {
        "ROI_{struct}".format(struct=structure).upper(): structure_index
        for structure_index, structure in enumerate(structure_names)
    }

    ROI_INV_DICT = {
        structure_index: "ROI_{struct}".format(struct=structure).upper()
        for structure_index, structure in enumerate(structure_names)
    }

    cont_map = [
        [
            [] for __ in range(SIZE_OF_GRID[2])
        ] for __ in range(0, len(rt.ROIContourSequence))
    ]

    for index in range(len(cont_map)):

        if index in ROI_INV_DICT:
            file_path = (root_dir 
                + '{structure}'.format(structure=ROI_INV_DICT[index]) 
                + '_contour.txt')
            if isfile(file_path):
                with open(file_path, 'rb') as struct_file:
                    cont_map[index] = pload(struct_file)
            else:
                print(
                    "No {structure} contour".format(structure=ROI_INV_DICT[index])
                    )
        else:
            for index2 in range(SIZE_OF_GRID[2]):
                cont_map[index][index2] = False
        
    print("\r" + " Finished ROI # {0}".format(index))

    return cont_map

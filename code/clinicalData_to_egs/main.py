#!/usr/bin/env
# author: Joseph Lucero
# created on: 27 May 2018 15:39:30
# purpose: turning patient CT files and dicom inputs to egs_brachy compatible files

# import system type modules
from __future__ import division

from sys import argv, stdout
from os import listdir, walk, getcwd, makedirs
from os.path import join, isfile, isdir
from pickle import load, dump
from time import time
from datetime import date
from copy import deepcopy

# load self-made modules
from process import process_CT, process_RP, process_RS
from get_funcs import get_target_dir, get_ref_info_from_ref_slice
from file_functions import create_EGSPhant
from special_funcs import perform_metallic_artifact_reduction as MAR

def main(egsphant_create=True, scale_dose=False, run_metallic_reduction=False):

    """\
    Description: Main function of code, calls other function to construct \
    egsphant from 

    Inputs:
    :param egsphant_create: Switch to run egsphant creation
    :type egsphant_create: bool
    :param scale_dose: Switch to scale dose or not
    :type scale_dose: bool
    :param run_metallic_reduction: Switch to run metallic artifiact reduction or not
    :type run_metallic_reduction: bool

    Outputs:
    VOID\
    """

    root_dir = get_target_dir()
    cwd = getcwd()
    egsphant_name = str(input(
        "Please input the base name of the egsphant:>> "
    ))
    
    ordered_CT_lst, slice_array_lst, num_CT_files = process_CT(root_dir)
    final_seed_locations, DSF = process_RP(root_dir, scale_dose)
    reference_dict = get_ref_info_from_ref_slice(
        *ordered_CT_lst[:2], num_CT_files=num_CT_files
    )
    cont_map = process_RS(root_dir, reference_dict["SIZE OF GRID"])
    
    if run_metallic_reduction:
        slice_array_lst = MAR(
            final_seed_locations, slice_array_lst, reference_dict["BOUNDS"],
            reference_dict["VOXEL DIMS"], reference_dict["VOXEL CENTERS"], 
            reference_dict["INTERCEPT"], reference_dict["SIZE OF GRID"]
            )

    if egsphant_create:
        create_EGSPhant(
            slice_array_lst=slice_array_lst, 
            contour_map=cont_map,
            ref_intercept=reference_dict["INTERCEPT"], 
            SIZE_OF_GRID=reference_dict["SIZE OF GRID"], 
            bounds=reference_dict["BOUNDS"],
            path_to_calibration='ask',
            path_to_media='ask',
            path_to_egsphants=None,
            egsphant_name=egsphant_name
        )

    return 0

if __name__ == "__main__":
    main()

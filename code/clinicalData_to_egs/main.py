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

# load external module functions
from dicom import read_file
from numpy import interp

# load self-made modules
from process import process_CT, process_RP
from get_funcs import get_target_dir
from special_funcs import perform_metallic_artifact_reduction as MAR

def main(scale_dose=False):

    """
    Description:

    Inputs:

    Outputs:
    """

    target_dir = get_target_dir()
    
    full_file_path = join(target_dir, listdir(target_dir)[0]) + '/'
    dir_sort = listdir(target_dir)

    for current_dir in dir_sort:
        if 'CT' in current_dir:
            file_path_CT = full_file_path + current_dir + '/'
        elif 'RP' in current_dir:
            file_path_RT = target_dir + current_dir + '/'
        elif 'RT' in current_dir:
            file_path_CONTOUR = target_dir + current_dir + '/'

    ordered_CT_lst, slice_array_lst = process_CT(target_dir)

    print ordered_CT_lst, slice_array_lst

    
if __name__ == "__main__":
    main()

#!/usr/bin/env
# author: Joseph Lucero
# created on: 25 May 2018 22:13:43
# purpose: define functions that handles files

# import external functions
from os import getcwd
from os.path import exists
from glob import glob
from numpy import array, empty

from get_funcs import get_target_dir()

def load_phantom(patientID, label, p_width):
    """
    Description: 
    Function that takes in patient ID and label to return the 
    doses and errors associated with that patient in 1D arrays.

    Input:
    :param patientID:
    :type patientID:
    :param label:
    :type label:
    :param p_width:
    :type p_width:

    Output:
    :param dose_array:
    :type dose_array:
    :param error_array:
    :type error_array:
    """

    p_data = []
    comp_data = []
    tissue_names = empty(p_width, dtype=int)
    target_dir = get_target_dir()

    while True:
        phantom_filename = raw_input(
            "Please input path to target file: "
        )
        if not exists(target_dir):
            print "That was not a valid path. " 
            continue
        else:
            break
    
    with open(phantom_filename,'r+') as c_phantom_file:
        for line in c_phantom_file:
            p_data.append(line)
    
    num_of_tissues = int(p_data)
    for tissue_index in xrange(num_of_tissues):
        tissue_names.append(p_data[tissue_index].rstrip("\n"))
    
    for line in p_data:
        p_line = line.rstrip("\n")

        # Tissue assignment lines in EGSphants are a big string equal to size of
        # phantom width
        if len(p_line) == p_width:
            for index, digit in enumerate(p_line):
                comp_data[index] = int(digit) 

    return tissue_names, comp_data
            
def create_EGSPhant():

    """
    Description: 
    Script to create EGSPhant files from CT files given a specific 
    CT calibration curve.

    Inputs:
    :param :
    :type :

    Outputs:
    :param :
    :type :
    """

    print "\n\n" + "*"*50
    print "Attempting to create egsphant files..."
    print "*"*50

    
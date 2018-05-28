#!/usr/bin/env
# author: Joseph Lucero
# created on: 27 May 2018 15:39:30
# purpose: turning patient CT files and dicom inputs to egs_brachy compatible files

# import system type modules
from sys import argv, stdout
from os import listdir, walk
from os.path import join, isfile, isdir, makedirs
from pickle import load, dump
from time import time
from datetime import date
from copy import deepcopy

# load external module functions
from dicom import read_file
from numpy import interp

def get_defaults():

    """
    Description:
    Define some default settings that are needed for the code

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

    return (
        NAME_STRING, NUM_OF_HISTORIES, AG_CUTOFF_DEN, CALC_CUTOFF_DEN, 
        BONE_CUTOFF_DEN, LOW_ARTIFACT_CUTOFF_DEN, CARTILAGE_CUTOFF_DEN,
        YELLOW_MARROW_CUTOFF_DEN, RED_MARROW_CUTOFF_DEN
        )

def main():

    """
    Description:

    Inputs:

    Outputs:
    """

    while True:
        target_dir = raw_input(
            """
            Please input path to directory in which all relevant files are stored :>> 
            """
        )
        if not isdir(target_dir):
            print "That is not a valid path to a directory."
            continue
        else:
            break

    patient_numbers = raw_input(
        "Please enter 2-digit patient number(s). (Delimit by spaces) :>> "
    )
    start_time = time()

    half_cutoffs = []
    with open(target_dir + '/half_cutoffs.txt','rb') as cutoff_file:
        for line in cutoff_file:
            half_cutoffs.append(float(line))
    half_counter = 0

    if patient_numbers == "all":
        pts = listdir(target_dir + "/MC_PBSI_COMPLETE")
        patient_numbers = []

        for pt in pts:
            patient_numbers.append(pts.split("Pt ")[1])
    else:
        patient_numbers = patient_numbers.split()

    



    
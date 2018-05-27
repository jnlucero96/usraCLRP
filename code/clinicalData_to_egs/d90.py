#!/usr/bin/env
# author: Joseph Lucero
# created on: 25 May 2018 22:13:43
# purpose: 

# load external functions 
from os import listdir, walk
from os.path import join, isfile, exists
from time import time
from dicom import read_file
from glob import glob
from numpy import array
from matplotlib.pyplot import subplots
from pickle import load, dump

# load self-modules
from py3ddose import DoseFile as read_dose

def plot_D90():

    """

    """

def get_CT(target_dir, patient_number, label, working_contours, no_skip=True):

    """
    Description: 
    Function that gets patient CT information. Opens, sorts, and arrays the
    CT files.

    Inputs:
    :param target_dir:
    :type target_dir:
    :param patient_numbers:
    :type patient_numbers: 
    """

    ct_filename_list = []
    sort_dcm_list = []
    ordered_ct_list = []

    dir_path = target_dir + "/Pt_%s/" % patient_number
    patient_file_path = join(dir_path,listdir(dir_path)[0]) + "/"

    for current_dir in listdir(patient_file_path):
        if "_CT_" in current_dir:
            ct_file_path = patient_file_path + current_dir + "/"
        elif "_RTPLAN_" in current_dir:
            plan_file_path = patient_file_path + current_dir + "/"
        elif "_RTst_" in current_dir:
            contour_file_path = patient_file_path + current_dir + "/"
    
    print "Now importing and arranging CT data"
    print "=" * 50

    for __, __, filenames in walk(ct_file_path):  
        ct_filename_list.extend(filenames)
        break  # limit search to first directory
    
    for ct_file in ct_filename_list: 
        # for loop to remove .mim files
        if 'dcm' not in ct_file:
            ct_filename_list.remove(ct_file)
            break  #JNL: not sure why this is here. Going to comment it out 4 now.
        else:
            path_to_read = ct_file_path + ct_file
            sort_dcm_list.append(
                (
                    read_file(path_to_read).SliceLocation,
                    path_to_read
                )
            )
    sort_dcm_list.sort()  # sort from inferior (-) to superior (+)
    
    for __, dicom_file in sort_dcm_list:
        ordered_ct_list.append(dicom_file)

    ref = read_file(ordered_ct_list[0]) # load reference file
    grid_size = (int(ref.Columns), int(ref.Rows), len(ordered_ct_list)) # order: (x,y,z)
    r_sog = (grid_size[2], grid_size[1], grid_size[0]) # order: (z,y,x) needed order for reshape

    try:
        dose_filename = glob(
            target_dir + '/3ddose/%s/*Pt_%s*' % (label, patient_number)
            )
        mc_dose_array = read_dose(dose_filename,load_uncertainty=False).dose
    except:
        no_skip = True

    if no_skip:

        for i, working_contour in enumerate(working_contours):

            cont_dose_array = []
            cont_dose_array_full = []

            if isfile()







def main():

    while True:
        target_dir = raw_input(
            "Please input top directory where relevant files are located :>> "
            )
        if not exists(target_dir):
            print "Not a valid path."
            continue
        else:
            break

    patient_numbers = raw_input(
        """
        Please enter 2-digit patient number(s). (Delimit entries by spaces) :>> 
        """
        )
    naming_string = raw_input(
        """
        Please input name of 3DDose files to use :>> 
        """
    )
    contours = raw_input(
        """
        Which contour(s)? (Delimit entries by spaces): lung, heart, skin_margin, 
        ptv1, ptv05, ctv, skin_skin, skin_surface, ribs, chest_wall, breast,
        body :>> 
        """
    ).split()

    start_time = time()

    if patient_numbers == 'all':
        patient_numbers = []
        write_trigger = True
        pts = listdir(target_dir + '/MC_PBSI_COMPLETE')

        for directory in listdir(target_dir):
            patient_numbers.append(directory.split('Pt ')[1])
    else:
        patient_numbers = patient_numbers.split()

    working_contours = working_contours.split()

    d90_list = [
        [] for i in xrange(len(working_contours))
    ]

    for patient_index in patient_numbers:
        print "\n*****************************************"
        print "******* Now working with Patient %s *******" % (patient_index)
        print "*****************************************\n"
    
    

    



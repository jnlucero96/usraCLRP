#!/usr/bin/env
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

    ###########################################################################
    ########
    ######## Opening, sorting, and arraying the CT files
    ########
    ###########################################################################

    print '\n\n' + '=' * 50
    print "Now importing and arranging CT data"
    print '=' * 50

    # extract the CT filenames in the target directory
    filename_CT_lst = glob(join(root_dir, 'CT*.dcm'))

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

    ###########################################################################
    ########
    ######## Open plan dicom files and process them
    ########
    ###########################################################################

    # extract the filenames in the RT directory
    filename_RP_lst = glob(join(root_dir, 'RP*.dcm')); answer = 0 

    num_RP_files = len(filename_RP_lst)

    # there should only be a unique RT (RP) plan file for a CT directory.
    # this covers the cases when user has more than one in the directory 
    if num_RP_files > 1:
        print "There is more than one RT plan in this directory:"
        for index, file_name in enumerate(filename_RP_lst):
            print str(index) + ": " + file_name
        answer = int(raw_input(
            "Which of the {0} files do you wish to use as the plan file? Please input in the index [int]:>> ".format(index + 1)
        ))

    current_plan_file = filename_RP_lst[answer]

    rp = dcmread(current_plan_file)
    print "Now finding treatement plan info. for Patient {0}.".format(
        rp.PatientID
    )

    try:
        air_kerma = float(rp.SourceSequence[0].ReferenceAirKermaRate)
        print 'Air Kerma =', air_kerma
    except (NameError, AttributeError):
        print "No Air Kerma found for Patient {0}.".format(rp.PatientID)
        scale_dose = False  # without air kerma cannot scale dose

    # define dose scaling factor (DSF)
    if scale_dose:
        # what is this hard coded number? for what source?
        # from Stephen this is for Pd-103
        DSF = 91558000000000 * air_kerma
    else:
        DSF = 1.0
        print "Dose remains unscaled."

    try:
        # get final list of seed locations
        final_seed_locations = [
            [
                float(
                    seed_info.ChannelSequence[0].BrachyControlPointSequence[1].ControlPoint3DPosition[i])
                for i in xrange(0, 3)
            ] for __, seed_info in enumerate(rp.ApplicationSetupSequence)
        ]
        print "Seed locations found."
    except:
        print "No seed locations."

    return final_seed_locations, DSF

def process_RS(root_dir, SIZE_OF_GRID):

    # glob the contour files in the directory
    contour_files = glob(join(root_dir, 'RS*.dcm')); answer = 0

    num_RS_files = len(contour_files)
    
    # there should only be a unique contour (RS) file for a CT directory.
    # this covers the cases when user has more than one in the directory
    if num_RS_files > 1:
        print "There is more than one RS plan in this directory:"
        for index, file_name in enumerate(contour_files):
            print str(index) + ": " + file_name
        answer = int(raw_input(
            "Which of the {0} files do you wish to use as the plan file? Please input in the index [int]:>> ".format(
                index + 1)
        ))

    current_contour_file = contour_files[answer]

    rt = dcmread(current_contour_file)

    structure_names = [
        name.ROIName for name in rt.StructureSetROISequence
    ]

    ROI_LUNG = -1
    ROI_HEART = -1 
    ROI_SKIN = -1
    ROI_SKIN_MARGIN = -1
    ROI_SKIN_SURFACE = -1
    ROI_PTV = -1
    ROI_CTV = -1
    ROI_BREAST = -1
    ROI_RIBS = -1
    ROI_CHEST_WALL = -1
    ROI_BODY = -1

    for relevant_index, roi_gen in enumerate(structure_names):
        ROI = roi_gen.lower()
        if ROI in ('ipsilateral lung', 'ipsilateral lungs'):
            ROI_LUNG = relevant_index
        elif ROI == 'heart':
            ROI_HEART = relevant_index
        elif ROI == 'skin margin':
            ROI_SKIN_MARGIN = relevant_index
        elif ROI == 'breast':
            ROI_BREAST = relevant_index
        elif ROI in ('ctv', 'final ctv df'):
            ROI_CTV = relevant_index
        elif ROI == 'skin':
            ROI_SKIN = relevant_index
        elif ROI == 'skin surface':
            ROI_SKIN_SURFACE = relevant_index
        elif ROI == 'ribs':
            ROI_RIBS = relevant_index
        elif ROI in ('chestwall', 'chest wall'):
            ROI_CHEST_WALL = relevant_index
        elif ROI == 'body':
            ROI_BODY = relevant_index
        else:
            pass

    if ROI_PTV == -1:  # these cases are misnamed. Catch them.
        for relevant_index, roi_gen in enumerate(structure_names):
            ROI = roi_gen.lower()
            if ROI in ('etv 1.0', 'ctv'):
                ROI_PTV = relevant_index
                if ROI == 'etv 1.0':
                    print '*'*10 + 'WARNING: Using ETV Structure' + '*'*10
                elif ROI == 'ctv':
                    print '*'*10 + 'WARNING: Using CTV Structure' + '*'*10
                else:
                    pass
            else: 
                pass

    if ROI_LUNG == -1:
        print "*"*10 + "WARNING: Could not find Lung Structure" + "*"*10
    if ROI_HEART == -1:
        print "*"*10 + "WARNING: Could not find Heart Structure" + "*"*10
    if ROI_SKIN_MARGIN == -1:
        print "*"*10 + "WARNING: Could not find Skin Structure" + "*"*10
    if ROI_PTV == -1:
        print "*"*10 + "WARNING: Could not find PTV Structure" + "*"*10
    if ROI_CTV == -1:
        print "*"*10 + "WARNING: Could not find CTV Structure" + "*"*10
    if ROI_SKIN == -1:
        print "*"*10 + "WARNING: Could not find Skin Structure" + "*"*10
    if ROI_SKIN_SURFACE == -1:
        print "*"*10 + "WARNING: Could not find Skin Surface Structure" + "*"*10

    cont_map = [
        [
            [] for __ in xrange(SIZE_OF_GRID[2])
        ] for __ in xrange(0, len(rt.ROIContourSequence))
    ]

    for index in xrange(len(cont_map)):

        if index == ROI_HEART:
            if isfile(root_dir + 'heart_contour.txt'):
                with open(root_dir + 'heart_contour.txt', 'rb') as heart_file:
                    cont_map[index] = pload(heart_file)
            else:
                print "No heart contour"
        elif index == ROI_LUNG:
            if isfile(root_dir + 'lung_contour.txt'):
                with open(root_dir + 'lung_contour.txt', 'rb') as lung_file:
                    cont_map[index] = pload(lung_file)
            else:
                print "No lung contour"
        elif index == ROI_BREAST:
            if isfile(root_dir + 'breast_contour.txt'):
                with open(
                    root_dir + 'breast_contour.txt', 'rb'
                    ) as breast_file:
                    cont_map[index] = pload(breast_file)
            else:
                print "No breast contour"
        elif index == ROI_SKIN:
            if isfile(root_dir + 'skin_skin_contour.txt'):
                with open(
                    root_dir + 'skin_skin_contour.txt', 'rb'
                    ) as skin_skin_file:
                    cont_map[index] = pload(skin_skin_file)
            else:
                print "No Actual Skin Contour"
        elif index == ROI_RIBS:
            if isfile(root_dir + 'ribs_contour.txt'):
                with open(root_dir + 'ribs_contour.txt', 'rb') as ribs_file:
                    cont_map[index] = pload(ribs_file)
        elif index == ROI_CHEST_WALL:
            if isfile(root_dir + 'chest_wall_contour.txt'):
                with open(
                    root_dir + 'chest_wall_contour.txt', 'rb'
                    ) as chest_wall_file:
                    cont_map[index] = pload(chest_wall_file)
            else:
                print "No chest wall contour"
        else:
            for index2 in xrange(SIZE_OF_GRID[2]):
                cont_map[index][index2] = False
        
    print "\r" + " Finished ROI # {0}".format(index)

    return cont_map

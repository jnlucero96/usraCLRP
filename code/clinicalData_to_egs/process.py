#!/usr/bin/env
# author: Joseph Lucero
# created on: 28 June 2018 17:43:28
# purpose: functions that process files and extract information

# import python modules
from __future__ import division

from os import walk
from pickle import load as pload
from os.path import join, isfile
from dicom import read_file


def process_CT(file_path):

    ###########################################################################
    ########
    ######## Opening, sorting, and arraying the CT files
    ########
    ###########################################################################

    print "Now importing and arranging CT data"
    print '=' * 50

    # extract the CT filenames in the target directory
    (__, __, filename_CT_lst) = walk(file_path).next()

    sort_dcm = [
        (
            read_file(file_path + CT_file).Slicelocation,
            file_path + CT_file
        ) for CT_file in filename_CT_lst
    ].sort()  # sort the dcm files from inferior (-) to superior (+)

    ordered_CT_lst = [
        item[1] for item in sort_dcm
    ]  # extract the ordered CT Dicom files

    slice_array_lst = [
        read_file(current_CT).pixel_array for current_CT in ordered_CT_lst
    ]  # slice array list indexing goes [z][y][x]

    return ordered_CT_lst, slice_array_lst


def process_RT(file_path, scale_dose):

    ###########################################################################
    ########
    ######## Open plan dicom files and process them
    ########
    ###########################################################################

    # extract the filenames in the RT directory
    root_dir, __, filename_lst = walk(file_path).next()

    for filename in filename_lst:
        if 'RTPLAN' and 'dcm' in join(root_dir, filename):
            current_plan_file = join(root_dir, filename)

    rp = read_file(current_plan_file)
    print "Now finding treatement plan info. for Patient {0}.".format(
        rp.PatientID
    )

    try:
        air_kerma = float(rp.SourceSequence[0].ReferenceAirKermaRate)
        print 'Air Kerma =', air_kerma
    except:
        print "No Air Kerma found for Patient {0}.".format(rp.PatientID)
        scale_dose = False  # without air kerma cannot scale dose

    # define dose scaling factor (DSF)
    if scale_dose:
        try:
            # what is this hard coded number? for what source?
            DSF = 91558000000000 * air_kerma
        except NameError:
            print "No air kerma info --> Cannot scale dose."
    else:
        DSF = 1
        print "Not scaling dose according to user preference."

    try:
        # get final list of seed locations
        final_seed_locations = [
            [
                float(seed_info.ChannelSequence[0].
                      BrachyControlPointSequence[1].ControlPoint3DPosition[i])
                for i in xrange(0, 3)
            ] for __, seed_info in enumerate(rp.ApplicationSetupSequence)
        ]
        print "Seed locations found."
    except:
        print "No seed locations found for Patient {0}.".format(rp.PatientID)

    return final_seed_locations, DSF

def process_CONTOUR(file_path, SIZE_OF_GRID):

    root_dir, __, filenames = walk(file_path).next()
    
    for file_name in filenames:
        if 'RTst' and 'dcm' in join(root_dir, file_name):
            current_contour_file = join(root_dir, file_name)

    rt = read_file(current_contour_file)

    structure_names = [
        name.ROI_name for name in rt.StructureSetROISequeunce
    ]

    print "Structures Found: "
    print structure_names

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
            if isfile(file_path + 'heart_contour.txt'):
                with open(file_path + 'heart_contour.txt', 'rb') as heart_file:
                    cont_map[index] = pload(heart_file)
            else:
                print "No heart contour"
        elif index == ROI_LUNG:
            if isfile(file_path + 'lung_contour.txt'):
                with open(file_path + 'lung_contour.txt', 'rb') as lung_file:
                    cont_map[index] = pload(lung_file)
            else:
                print "No lung contour"
        elif index == ROI_BREAST:
            if isfile(file_path + 'breast_contour.txt'):
                with open(
                    file_path + 'breast_contour.txt', 'rb'
                    ) as breast_file:
                    cont_map[index] = pload(breast_file)
            else:
                print "No breast contour"
        elif index == ROI_SKIN:
            if isfile(file_path + 'skin_skin_contour.txt'):
                with open(
                    file_path + 'skin_skin_contour.txt', 'rb'
                    ) as skin_skin_file:
                    cont_map[index] = pload(skin_skin_file)
            else:
                print "No Actual Skin Contour"
        elif index == ROI_RIBS:
            if isfile(file_path + 'ribs_contour.txt'):
                with open(file_path + 'ribs_contour.txt', 'rb') as ribs_file:
                    cont_map[index] = pload(ribs_file)
        elif index == ROI_CHEST_WALL:
            if isfile(file_path + 'chest_wall_contour.txt'):
                with open(
                    file_path + 'chest_wall_contour.txt', 'rb'
                    ) as chest_wall_file:
                    cont_map[index] = pload(chest_wall_file)
            else:
                print "No chest wall contour"
        else:
            for index2 in xrange(SIZE_OF_GRID[2]):
                cont_map[index][index2] = False
        
    print "\r" + " Finished ROI # %s" % str(index)
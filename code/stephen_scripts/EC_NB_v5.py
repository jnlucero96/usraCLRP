#!/usr/bin/env

import sys
import dicom
import numpy as np
import os
from os import walk
from matplotlib.path import Path
import pickle
import time
import datetime
import copy

## INDIVIDUALIZATION IS CURRENTLY TURNED OFF. FIX DENSITIES BEFORE TURNING ON ## THink they've been fixed? Check this

########################
Name_String = 'Af_Final'
Number_of_Histories = '1e9'
AG_Cutoff_Den = 0.985
Calc_Cutoff_Den = 1.16
Bone_Cutoff_Den = 1.16
Low_Artifact_Cutoff_Den = 0.80
Cartilage_Cutoff_Den = 1.066
Yellow_Marrow_Cutoff_Den = 1.005
Red_Marrow_Cutoff_Den = 1.055

Indiv = False # Turns On/Off Individualization
########################

PatientNumbers = raw_input('Please enter 2-Digit Patient Number(s) (Delimited by Spaces) #: ')
start_time = time.time()

half_cutoffs = []
cutoff_file = open('./half_cutoffs.txt', 'rb')
for line in cutoff_file:
    half_cutoffs.append(float(line))
half_counter = 0
#print(half_cutoffs)

# Comes from Breast_Composition_Individual_Cut_Compare.py

# Getting individual Gland/Adipose Boundary

# density_path = 'density.txt'
# densities = open(density_path, 'r')
# adipose_peaks = [0]
# gland_peaks = [0]
#
# for lines in densities:
#     if 'Adipose Gaussian Peak' in lines:
#         adipose_peaks.append(float(lines[38:45]))
#     if 'Gland' in lines:
#         gland_peaks.append(float(lines[27:33]))

#print('Adipose Peaks are at:')
#print(adipose_peaks)
#print('Gland Peaks are at:')
#print(gland_peaks)



if PatientNumbers == 'all':
    pts = os.listdir('./MC PBSI Complete')
    PatientNumbers = []

    for i in pts:
        PatientNumbers.append(i.split('Pt ')[1])

else:
    PatientNumbers = PatientNumbers.split(' ')

for PatientNumber in PatientNumbers:
    print('\n*************************************')
    print('***  Now working with Patient %s  ***' % PatientNumber)
    print('*************************************\n')

    FirstFilePath = './MC PBSI Complete/Pt %s/' % PatientNumber
    FilePath = os.path.join(FirstFilePath, os.listdir(FirstFilePath)[0]) + '/'
    dir_sort = os.listdir(FilePath)
    for current_dir in dir_sort:
        if '_CT_' in current_dir:
            CTFilePath = FilePath + current_dir + '/'
        elif '_RTPLAN_' in current_dir:
            PlanFilePath = FilePath + current_dir + '/'
        elif '_RTst_' in current_dir:
            ContourFilePath = FilePath + current_dir + '/'

    ############################################################################
    # Opening and Sorting and Arraying the CT Files
    ############################################################################

    print('Importing and Arranging CT Data')
    print('-' * 50)
    ctfilenamelist = []
    for dirname, dirnames, filenames in walk(CTFilePath):
        ctfilenamelist.extend(filenames)
        break

    sort_dc = []
    for ct_file in ctfilenamelist:  # Removes .mim files from list
        if 'dcm' not in ct_file:
            ctfilenamelist.remove(ct_file)
            break  # I think this break needs to go if there's more than just .mim files

        sort_dc.append((dicom.read_file(CTFilePath + ct_file).SliceLocation, CTFilePath + ct_file))

    sort_dc.sort()  # Sorts the DCM Files from Inferior (-) to Superior (+)
    OrderedCTList = []
    for slice_location, DIFile in sort_dc:
        OrderedCTList.append(DIFile)

    # Now have ordered list of CT Dicom files for specified patient.
    # In the SliceListArray it goes [z][y][x]

    SliceListArray = []

    for current_ct in OrderedCTList:
        # print('Starting a New Slice')
        ct = dicom.read_file(current_ct)
        SliceListArray.append(ct.pixel_array)

    #######################################################################
    # Open Plan Dicom Files
    #######################################################################

    for root, dirnames, filenames in os.walk(PlanFilePath):
        for filename in filenames:
            if 'RTPLAN' in os.path.join(root, filename) and 'dcm' in os.path.join(root, filename):
                current_plan_file = os.path.join(root, filename)

    rp = dicom.read_file(current_plan_file)
    print('Finding Treatment Plan Information for Patient %s' % str(rp.PatientID))

    # Extract the Air Kerma Strength of the Seeds.
    try:
        AirKerma = float(rp.SourceSequence[0].ReferenceAirKermaRate)
        # print(AirKerma)
        print('Air Kerma = %s' % str(AirKerma))
    except:
        print('No Air Kerma found for Patient # %d' % PatientNumber)

    DSF = 91558000000000 * AirKerma  # Dose Scaling Factor = (tau(Pd-103)/100) * (AirKerma / EGS Brachy Theraseed Air Kerma per History)  -- 100 is for cGy -> Gy

    # Extract the Seed Locations // Make Full Seed Transformation File Creation Here with a While() statement to turn on/off
    try:
        SeedLocations = []
        for seednumber, seedinfo in enumerate(rp.ApplicationSetupSequence):
            SeedLocation = []

            for i in range(0, 3):
                SeedLocation.append(
                    float(seedinfo.ChannelSequence[0].BrachyControlPointSequence[1].ControlPoint3DPosition[i]))

            # print(seednumber+1,SeedLocation)
            SeedLocations.append(SeedLocation)  # Final List of Seed Locations

        print('Seed Locations Found')
    except:
        print('No Seed Locations Found for Patient # %d') % int(rp.PatientID)

        # Get Essential Voxel and Slice Information from Reference Slice (First Slice)/

    Ref = dicom.read_file(OrderedCTList[0])
    Ref2 = dicom.read_file(OrderedCTList[1])
    Slice_Thickness = Ref2.SliceLocation - Ref.SliceLocation
    Start_of_Grid = (float(Ref.ImagePositionPatient[0]), float(Ref.ImagePositionPatient[1]), float(Ref.SliceLocation))
    Size_of_Grid = (int(Ref.Columns), int(Ref.Rows), len(OrderedCTList))
    Voxel_Size = (float(Ref.PixelSpacing[0]), float(Ref.PixelSpacing[1]), Slice_Thickness)
    Intercept = float(Ref.RescaleIntercept)
    Scaling = float(Ref.RescaleSlope)
    print(Voxel_Size)

    # Set Bounds for Voxels

    xbounds = []
    xbounds.append((Start_of_Grid[0] - 0.5 * Voxel_Size[0]) / 10)
    for i in range(Size_of_Grid[0]):
        xbounds.append(xbounds[i] + Voxel_Size[0] / 10)

    ybounds = []
    ybounds.append((Start_of_Grid[1] - 0.5 * Voxel_Size[1]) / 10)
    for i in range(Size_of_Grid[1]):
        ybounds.append(ybounds[i] + Voxel_Size[1] / 10)

    zbounds = []
    zbounds.append((Start_of_Grid[2] - 0.5 * Voxel_Size[2]) / 10)
    for i in range(Size_of_Grid[2]):
        zbounds.append(zbounds[i] + Voxel_Size[2] / 10)

    voxel_centers = []
    for i in range(len(ybounds) - 1):
        for j in range(len(xbounds) - 1):
            voxel_centers.append(((xbounds[j] + 0.5 * Voxel_Size[0] / 10), (ybounds[i] + 0.5 * Voxel_Size[1] / 10)))

    ###################################################################
    ## Metallic Artifact Reduction
    ###################################################################

    print('\nAttempting to preform Metallic Artifact Reduction')
    print('-' * 50)

    Replace = -90.2
    Cutoff = 200  # In BIG STR I lowered cutoff 200 -> 150 - > 100, increased MaxRadius to 0.5 -> 1.0, and Went from +/-2 slides to +/- 4
    CylMaxRadius = 0.5
    MAR_counter = 0

    for Seed in SeedLocations:

        sliceheight = int((Seed[2] - 10 * zbounds[0]) // Voxel_Size[2])  # Gives the index for the slice
        # index starts at 0
        for index, center in enumerate(voxel_centers):
            xdis = abs(center[0] - Seed[0] / 10) - 0.05 * Voxel_Size[0]
            ydis = abs(center[1] - Seed[1] / 10) - 0.05 * Voxel_Size[1]

            if (xdis ** 2 + ydis ** 2) < CylMaxRadius ** 2:  # This voxel is inside the checked cylinder
                xpos = index % Size_of_Grid[0]
                ypos = index // Size_of_Grid[0]

                if SliceListArray[sliceheight][ypos][xpos] + Intercept > Cutoff:
                    MAR_counter += 1
                    SliceListArray[sliceheight][ypos][xpos] = Replace - Intercept

                if SliceListArray[sliceheight + 1][ypos][xpos] + Intercept > Cutoff:
                    MAR_counter += 1
                    SliceListArray[sliceheight + 1][ypos][xpos] = Replace - Intercept

                if SliceListArray[sliceheight + 2][ypos][xpos] + Intercept > Cutoff:
                    MAR_counter += 1
                    SliceListArray[sliceheight + 2][ypos][xpos] = Replace - Intercept

                if SliceListArray[sliceheight - 1][ypos][xpos] + Intercept > Cutoff:
                    MAR_counter += 1
                    SliceListArray[sliceheight - 1][ypos][xpos] = Replace - Intercept

                if SliceListArray[sliceheight - 2][ypos][xpos] + Intercept > Cutoff:
                    MAR_counter += 1
                    SliceListArray[sliceheight - 2][ypos][xpos] = Replace - Intercept

    print('Metallic Artifact Reduction completed. %s high density voxels replaced' % str(MAR_counter))
    ###################################################################
    # Get Contour/Structure Data
    ###################################################################

    # Finding correct Dicom plan file
    for root, dirnames, filenames in os.walk(ContourFilePath):
        for filename in filenames:
            if 'RTst' in os.path.join(root, filename) and 'dcm' in os.path.join(root, filename):
                current_contour_file = os.path.join(root, filename)

    rt = dicom.read_file(current_contour_file)

    print('\nFinding Contour Information for Patient %s' % str(rt.PatientID))
    print('-' * 50)

    # Choosing which Contour to Plot

    structure_names = [name.ROIName for name in rt.StructureSetROISequence]
    print('Structures Found: ')
    print(structure_names)

    roi_lung = -1
    roi_heart = -1
    roi_skin_margin = -1
    roi_ptv = -1
    roi_ctv = -1

    for ROI in structure_names:
        if ROI.lower() == 'ipsilateral lung' or ROI.lower() == 'ipsilateral lungs':
            roi_lung = structure_names.index(ROI)
        elif ROI.lower() == 'heart':
            roi_heart = structure_names.index(ROI)
        elif ROI.lower() == 'skin margin':
            roi_skin_margin = structure_names.index(ROI)
        elif ROI.lower() == 'ptv 1.0' or ROI.lower() == 'etv 1.0':
            roi_ptv1 = structure_names.index(ROI)
        elif ROI.lower() == 'ptv 0.5' or ROI.lower() == 'etv 0.5':
            roi_ptv05 = structure_names.index(ROI)
        elif ROI.lower() == 'breast':
            roi_breast = structure_names.index(ROI)
        elif ROI.lower() == 'ctv' or ROI.lower() == 'final ctv df':
            roi_ctv = structure_names.index(ROI)
        elif ROI.lower() == 'skin':
            roi_skin = structure_names.index(ROI)
        elif ROI.lower() == 'skin surface':
            roi_skin_surface = structure_names.index(ROI)
        elif ROI.lower() == 'ribs':
            roi_ribs = structure_names.index(ROI)
        elif ROI.lower() == 'chestwall' or ROI.lower() == 'chest wall':
            roi_chest_wall = structure_names.index(ROI)
        elif ROI.lower() == 'body':
            roi_body = structure_names.index(ROI)

    # A couple cases are misnamed, this catches them.
    if roi_ptv == -1:
        for ROI in structure_names:
            if ROI.lower() == 'etv 1.0':
                roi_ptv = structure_names.index(ROI)
                print('**********Warning: Using ETV Structure**********')
    if roi_ptv == -1:
        for ROI in structure_names:
            if ROI.lower() == 'ctv':
                roi_ptv = structure_names.index(ROI)
                print('**********Warning: Using CTV Structure**********')

    if roi_lung == -1:
        print('***********Warning: Could not find Lung Structure************')
    if roi_heart == -1:
        print('***********Warning: Could not find Heart Structure************')
    if roi_skin_margin == -1:
        print('***********Warning: Could not find Skin Structure************')
    if roi_ptv == -1:
        print('***********Warning: Could not find PTV Structure************')
    if roi_ctv == -1:
        print('***********Warning: Could not find CTV Structure************')
    if roi_skin == -1:
        print('***********Warning: Could not find Skin Structure************')
    if roi_skin_surface == -1:
        print('***********Warning: Could not find Skin Surface Structure************')


    # Finding contour map for voxel assignment
    print('Finding Contour Masks for Patient %s' % str(rt.PatientID))

    contmap = [[[] for z in range(Size_of_Grid[2])] for ROI in range(0, len(rt.ROIContourSequence))]

    for i in range(len(contmap)):

        if i == roi_heart:
            if os.path.isfile(FilePath + 'heart_contour.txt'):
                contmap[i] = pickle.load(open(FilePath + 'heart_contour.txt', 'rb'))

            else:
                print('No Heart Contour')

        elif i == roi_lung:
            if os.path.isfile(FilePath + 'lung_contour.txt'):
                contmap[i] = pickle.load(open(FilePath + 'lung_contour.txt', 'rb'))

            else:
                print('No Lung Contour')

        elif i == roi_breast:
            if os.path.isfile(FilePath + 'breast_contour.txt'):
                contmap[i] = pickle.load(open(FilePath + 'breast_contour.txt', 'rb'))

            else:
                print('No Breast Contour')

        elif i == roi_skin:
            if os.path.isfile(FilePath + 'skin_skin_contour.txt'):
                contmap[i] = pickle.load(open(FilePath + 'skin_skin_contour.txt', 'rb'))

            else:
                print('No Actual Skin Contour')

        elif i == roi_ribs:
            if os.path.isfile(FilePath + 'ribs_contour.txt'):
                contmap[i] = pickle.load(open(FilePath + 'ribs_contour.txt', 'rb'))

            else:
                print('No Ribs Contour')

        elif i == roi_chest_wall:
            if os.path.isfile(FilePath + 'chest_wall_contour.txt'):
                contmap[i] = pickle.load(open(FilePath + 'chest_wall_contour.txt', 'rb'))

            else:
                print('No Chest Wall Contour')

        else:
            for j in range(Size_of_Grid[2]):
                contmap[i][j] = False

        sys.stdout.write("\r" + ' Finished ROI # %s ' % str(i))
        sys.stdout.flush()


    #######################################################################################################################
    ## EGSPhant Creation / Tissue Assignment and Voxel Density
    #######################################################################################################################

    Final_Skin = copy.deepcopy(contmap[roi_skin])

    print('\n\nAttempting to Create Egsphant Files')
    print('-' * 50)

    SD_x = [-1000, 0.00, 61.9, 1000, 2000, 3000, 3100, 5000, 10000, 20000, 25000]  # Mass Density Curve - x
    SD_y = [0.001, 1.008, 1.073, 1.667, 2.3, 2.933, 2.999, 2.999, 7.365, 10, 10]  # - y

    if Indiv:
        AG_Cutoff_Den = half_cutoffs[half_counter]
        print('Individual AG Density Cutoff Used is # %d: %f g / cm-3' % (half_counter, AG_Cutoff_Den))
        half_counter += 1

    AG_Cutoff_HU = np.interp(AG_Cutoff_Den, SD_y, SD_x)
    print('A_G_Cutoff_HU = ' + str(AG_Cutoff_HU))

    Calc_Cutoff_HU = np.interp(Calc_Cutoff_Den, SD_y, SD_x)
    print('Calcification Cutoff HU = ' + str(Calc_Cutoff_HU))

    Bone_Cutoff_HU = np.interp(Bone_Cutoff_Den, SD_y, SD_x)
    print('Bone Cutoff HU = ' + str(Bone_Cutoff_HU))

    Low_Artifact_Cutoff_HU = np.interp(Low_Artifact_Cutoff_Den, SD_y, SD_x)
    print('Low Artifact Cutoff HU = ' + str(Low_Artifact_Cutoff_HU))

    Cartilage_Cutoff_HU = np.interp(Cartilage_Cutoff_Den, SD_y, SD_x)
    print('Low Artifact Cutoff HU = ' + str(Cartilage_Cutoff_HU))

    Red_Marrow_Cutoff_HU = np.interp(Red_Marrow_Cutoff_Den, SD_y, SD_x)
    print('Low Artifact Cutoff HU = ' + str(Red_Marrow_Cutoff_HU))

    Yellow_Marrow_Cutoff_HU = np.interp(Yellow_Marrow_Cutoff_Den, SD_y, SD_x)
    print('Low Artifact Cutoff HU = ' + str(Yellow_Marrow_Cutoff_Den))



    Media = ['AIR_TG43','MUSCLE2_WW86','CARTILAGE_WW86', 'SKIN2_WW86', 'GLAND2_WW86', 'ADIPOSE2_WW86','HEART_BLOODFILLED_WW86', 'LUNG_BLOODFILLED_WW86', 'CORTICAL_BONE_WW86','YELLOW_MARROW_WW86','RED_MARROW_WW86','CALCIFICATION_ICRU46']

    # Opening File
    phantpath = inputpath = './EGSphants/%s/Br_%s_Pt_%s.egsphant' % (Name_String, Name_String, str(PatientNumber))
    if not os.path.isdir('./EGSphants/' + Name_String):
        print('Making Egsphant Folder:' + './EGSphants/' + Name_String)
        os.makedirs('./EGSphants/' + Name_String)
    egsphant = open(phantpath, 'w')

    # Header Data
    egsphant.write(str(len(Media)) + '\n')
    for tissue in Media:
        egsphant.write(tissue + '\n')  # Write Tissues

    egsphant.write('  0.25')
    for i in range(len(Media) - 1):
        egsphant.write(' ' * 7 + '0.25')
    egsphant.write('\n')

    for Length in Size_of_Grid:
        egsphant.write(' ' * 2 + str(Length))  # Mehan has extra thing here, check if egsphant breaks

    # Write Boundaries for Voxels
    newline_counter = 0
    print('Writing Voxel Boundaries')
    for x in xbounds:
        if newline_counter % 3 == 0:
            egsphant.write('\n' + '   ')
        egsphant.write("{:.5f}".format(x))
        if newline_counter % 3 != 2:
            egsphant.write(' ' * 10)
        elif newline_counter % 3 == 2:
            egsphant.write(' ' * 7)
        newline_counter += 1
    #
    for y in ybounds:
        if newline_counter % 3 == 0:
            egsphant.write('\n' + '   ')
        egsphant.write("{:.5f}".format(y))
        if newline_counter % 3 != 2:
            egsphant.write(' ' * 10)
        elif newline_counter % 3 == 2:
            egsphant.write(' ' * 7)
        newline_counter += 1

    # print(zbounds)
    for z in zbounds:
        if newline_counter % 3 == 0:
            egsphant.write('\n' + '   ')
        egsphant.write("{:.5f}".format(z))
        if newline_counter % 3 != 2:
            egsphant.write(' ' * 10)
        elif newline_counter % 3 == 2:
            egsphant.write(' ' * 7)
        newline_counter += 1
    egsphant.write('\n')
    print('.....finished')

    # Writing Tissue Assignments

    bone_count = 0
    gland_count = 0
    adipose_count = 0
    soft_tissue_count = 0
    calc_count = 0
    skin_count = 0
    air_in_skin = 0
    r_marrow_count = 0
    y_marrow_count = 0
    cart_count = 0
    muscle_count = 0
    adi_replace = 0

    print('Writing Tissue Assignments')
    for z in range(Size_of_Grid[2]):
        for y in range(Size_of_Grid[1]):
            for x in range(Size_of_Grid[0]):

                # Tissue Assignment Scheme is Here
                ct_value = SliceListArray[z][y][x] + Intercept

                ## ADD RIBS and CHEST WALL CONTOURS TO AVOID SKIN ASSIGNMENT!?!

                if contmap[roi_lung][z][x + y * Size_of_Grid[0]]:
                    if ct_value < -400:
                        egsphant.write(str(Media.index('AIR_TG43') + 1))
                    else:
                        egsphant.write(str(Media.index('LUNG_BLOODFILLED_WW86') + 1))

                elif contmap[roi_heart][z][x + y * Size_of_Grid[0]]:
                    egsphant.write(str(Media.index('HEART_BLOODFILLED_WW86') + 1))

                elif contmap[roi_ribs][z][x + y * Size_of_Grid[0]]:
                    if ct_value > Red_Marrow_Cutoff_HU:
                        egsphant.write(str(Media.index('CORTICAL_BONE_WW86') + 1))
                        bone_count += 1
                    elif ct_value < Yellow_Marrow_Cutoff_HU:
                        egsphant.write('A')
                        y_marrow_count += 1
                    else:
                        egsphant.write('B')
                        r_marrow_count += 1

                elif contmap[roi_chest_wall][z][x + y * Size_of_Grid[0]]:
                    if ct_value > Bone_Cutoff_HU:
                        egsphant.write(str(Media.index('CORTICAL_BONE_WW86') + 1))
                        bone_count += 1
                    elif ct_value > Cartilage_Cutoff_HU:
                        egsphant.write(str(Media.index('CARTILAGE_WW86') + 1))
                        cart_count += 1
                    elif ct_value > AG_Cutoff_HU:
                        egsphant.write(str(Media.index('MUSCLE2_WW86') + 1))
                        adipose_count += 1
                    else:
                        egsphant.write(str(Media.index('ADIPOSE2_WW86') + 1))
                        adipose_count += 1

                elif contmap[roi_skin][z][x + y * Size_of_Grid[0]]:
                    if ct_value < -400:
                        egsphant.write(str(Media.index('AIR_TG43') + 1))
                        air_in_skin += 1
                        Final_Skin[z][x + y * Size_of_Grid[0]] = False
                    else:
                        egsphant.write(str(Media.index('SKIN2_WW86') + 1))
                        skin_count += 1

                elif contmap[roi_breast][z][x + y * Size_of_Grid[0]]:
                    if ct_value > Calc_Cutoff_HU:
                        egsphant.write('C')
                        calc_count += 1
                    elif ct_value < AG_Cutoff_HU:
                        egsphant.write(str(Media.index('ADIPOSE2_WW86') + 1))
                        adipose_count += 1
                    else:
                        egsphant.write(str(Media.index('GLAND2_WW86') + 1))
                        gland_count += 1

                else:
                    if ct_value > Bone_Cutoff_HU:
                        egsphant.write(str(Media.index('CORTICAL_BONE_WW86') + 1))
                        bone_count += 1
                    elif ct_value < -400:
                        egsphant.write(str(Media.index('AIR_TG43') + 1))
                    elif ct_value > AG_Cutoff_HU:
                        egsphant.write(str(Media.index('MUSCLE2_WW86') + 1))
                        muscle_count += 1
                    else:
                        # This checks for adipose voxels lying on the outside of the skin contour. Due to imperfect contouring?
                        try:
                            if ((SliceListArray[z][y][x+1] + Intercept) < -400) or ((SliceListArray[z][y][x-1] + Intercept) < -400) or ((SliceListArray[z][y+1][x] + Intercept) < -400) or ((SliceListArray[z][y-1][x] + Intercept) < -400):
                                if contmap[roi_skin][z][x + 1 + y * Size_of_Grid[0]] or contmap[roi_skin][z][x - 1 + y * Size_of_Grid[0]] or contmap[roi_skin][z][x + (y+1) * Size_of_Grid[0]] or contmap[roi_skin][z][x + (y-1) * Size_of_Grid[0]]:
                                    egsphant.write(str(Media.index('SKIN2_WW86') + 1))
                                    Final_Skin[z][x + y * Size_of_Grid[0]] == True
                                    adi_replace += 1
                                else:
                                    egsphant.write(str(Media.index('ADIPOSE2_WW86') + 1))
                                    adipose_count += 1

                            else:
                                egsphant.write(str(Media.index('ADIPOSE2_WW86') + 1))
                                adipose_count += 1
                        except:
                            egsphant.write(str(Media.index('ADIPOSE2_WW86') + 1))
                            adipose_count += 1
                            print('The try statement worked!')

            egsphant.write('\n')
        egsphant.write('\n')

    pickle.dump(Final_Skin, open(FilePath + 'skin_final_contour.txt', 'wb'))
    print('.....finished')
    print('Outside Adipose Replaced = ' + str(adi_replace))


    # Writing Voxel Densities (Convert from Electron to Mass)

    print('Writing Voxel Density Assignments')
    low_artifact_counter = 0
    newline_counter = 0
    for z in range(Size_of_Grid[2]):
        sys.stdout.write("\r" + str(z) + '/' + str(Size_of_Grid[2]))
        sys.stdout.flush()
        for y in range(Size_of_Grid[1]):
            for x in range(Size_of_Grid[0]):
                ct_value = SliceListArray[z][y][x] + Intercept  # Should there be a function that limits the density inside the PTV? It's in MAR
                current_density = np.interp(ct_value, SD_x, SD_y)
                if contmap[roi_breast][z][x + y * Size_of_Grid[0]] and current_density < Low_Artifact_Cutoff_Den:
                    current_density = Low_Artifact_Cutoff_Den
                    low_artifact_counter += 1

                if newline_counter % 3 == 0:
                    egsphant.write('\n' + '  ')
                egsphant.write("{:.15f}".format(current_density))

                if newline_counter % 3 != 2:
                    egsphant.write(' ' * 9)
                elif newline_counter % 3 == 2:
                    egsphant.write(' ' * 7)
                newline_counter += 1

            egsphant.write('\n')
            newline_counter = 0
        egsphant.write('\n')

    print('Patient Information Output to: ' + './Patient_Information/%s/%s_Patient_%s_Info.txt' % (Name_String, Name_String, str(PatientNumber)))
    print('.....finished')

    egsphant.close()

    #######################################################################################################################
    ## EGS Input File Creation
    #######################################################################################################################

    print('\nAttempting to Create Input File...')
    print('-' * 50)

    # Opening File
    inputpath = './Inputs/%s/Br_%s_Pt_%s.egsinp' % (Name_String, Name_String, str(PatientNumber))
    if not os.path.exists('./Inputs/' + Name_String):
        os.makedirs('./Inputs/' + Name_String)
        print('Making Inputs Folder:' + '/Inputs/' + Name_String)

    egsinp = open(inputpath, 'wb')

    # Header Data
    egsinp.write('#' * 100 + '\n' + '#' * 5 + '\n' + '#' * 5)
    egsinp.write('  This is an automatically generated Breast Patient Input File created on: %s \n' % str(datetime.date.today()))
    egsinp.write('#####  Information for this file can be found at: ' + './Patient_Information/%s/%s_Patient_%s_Info.txt' % (Name_String, Name_String, str(PatientNumber)))
    egsinp.write('#####  File created by: ' + str(sys.argv[0]))

    egsinp.write('#' * 5 + '\n' + '#' * 100 + '\n\n')

    # Body
    egsinp.write(':start run control:\n')
    egsinp.write('	ncase = %s\n' % Number_of_Histories)
    egsinp.write('  nbatch = 4\n')
    egsinp.write('  nchunk = 1\n')
    egsinp.write('	#calculation = combine\n')
    egsinp.write('	geometry error limit = 1000000\n')
    egsinp.write(':stop run control:\n\n')
    egsinp.write('#' + '-' * 100 + '\n\n')

    egsinp.write(':start run mode:\n')
    egsinp.write('	run mode = normal\n')
    egsinp.write(':stop run mode:\n\n')
    egsinp.write('#' + '-' * 100 + '\n\n')

    egsinp.write(':start media definition:\n')
    egsinp.write('    AE = 1.512\n')
    egsinp.write('    UE = 2.012\n')
    egsinp.write('    AP = 0.001\n')
    egsinp.write('    UP = 1.500\n')

    egsinp.write('    material data file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/media/material.dat\n\n')
    egsinp.write('    :start HEART_BLOODFILLED_WW86:\n')
    egsinp.write('        bremsstrahlung correction = NRC\n')
    egsinp.write('        density correction file = heart_blood-filled_icru_1986\n')
    egsinp.write('    :stop HEART_BLOODFILLED_WW86:\n\n')


    egsinp.write(':stop media definition:\n\n')
    egsinp.write('#' + '-' * 100 + '\n\n')

    egsinp.write(':start geometry definition:\n')
    egsinp.write('source geometries = seed\n')
    egsinp.write('phantom geometries = phantom\n')
    egsinp.write('simulation geometry = phant_w_seeds\n\n')

    egsinp.write('	:start geometry:\n')
    egsinp.write('        name = phantom\n')
    egsinp.write('        library = egs_glib\n')
    egsinp.write('        type = egsphant\n')
    egsinp.write('        egsphant file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/geometry/phantoms/%s/Br_%s_Pt_%s.egsphant\n' % (Name_String, Name_String, str(PatientNumber)))
    egsinp.write('        density file = /data/data062/sdeering/EGSnrc_Sept2016/pegs4/data/brachy_xcom_1.5MeV.pegs4dat\n')
    egsinp.write('    :stop geometry:\n\n')

    egsinp.write('	:start geometry:\n')
    egsinp.write('		name = seed\n')
    egsinp.write('		library = egs_glib\n')
    egsinp.write('		include file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/geometry/sources/Pd103/TheraSeed200/TheraSeed200.geom\n')
    egsinp.write('	:stop geometry:\n')

    egsinp.write('	:start geometry:\n')
    egsinp.write('		name = phant_w_seeds\n')
    egsinp.write('		library = egs_autoenvelope\n')
    egsinp.write('		base geometry = phantom\n\n')
    egsinp.write('		:start inscribed geometry:\n')
    egsinp.write('			inscribed geometry name = seed\n')
    egsinp.write('			:start transformations:\n')
    egsinp.write('				include file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/geometry/transformations/Breast_Patients/Transformations/seeds_xy45_Pt_%s\n' % str(PatientNumber))
    egsinp.write('			:stop transformations:\n')
    egsinp.write('			:start volume correction:\n')
    egsinp.write('				correction type = correct\n')
    egsinp.write('				density of random points (cm^-3) = 1E8\n')
    egsinp.write('				include file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/geometry/sources/Pd103/TheraSeed200/boundary.shape\n')
    egsinp.write('			:stop volume correction:\n')
    egsinp.write('		:stop inscribed geometry:\n')
    egsinp.write('	:stop geometry:\n\n')
    egsinp.write(':stop geometry definition:\n\n')
    egsinp.write('#' + '-' * 100 + '\n\n')

    egsinp.write(':start volume correction:\n\n')
    egsinp.write('    :start source volume correction:\n')
    egsinp.write('        correction type = none\n')
    egsinp.write('        density of random points (cm^-3) = 1E8\n')
    egsinp.write('        include file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/geometry/sources/Pd103/TheraSeed200/boundary.shape\n')
    egsinp.write('    :stop source volume correction:\n\n')
    egsinp.write(':stop volume correction:\n\n')
    egsinp.write('#' + '-' * 100 + '\n\n')

    egsinp.write(':start source definition:\n\n')
    egsinp.write('    :start source:\n')
    egsinp.write('        name = Thera200\n')
    egsinp.write('        library = egs_isotropic_source\n')
    egsinp.write('        charge = 0\n')
    egsinp.write('        include file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/geometry/sources/Pd103/TheraSeed200/TheraSeed200.shape\n\n')
    egsinp.write('        :start spectrum:\n')
    egsinp.write('            type = tabulated spectrum\n')
    egsinp.write('            spectrum file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/spectra/Pd103_TG43.spectrum\n')
    egsinp.write('        :stop spectrum:\n')
    egsinp.write('    :stop source:\n\n')
    egsinp.write('    :start transformations:\n')
    egsinp.write('        include file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/geometry/transformations/Breast_Patients/Transformations/seeds_xy45_Pt_%s\n' % str(PatientNumber))
    egsinp.write('    :stop transformations:\n\n')
    egsinp.write('    simulation source = Thera200\n')
    egsinp.write(':stop source definition:\n\n')
    egsinp.write('#' + '-' * 100 + '\n\n')

    egsinp.write(':start scoring options:\n')
    egsinp.write('    score energy deposition = no\n')
    egsinp.write('    muen file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/muen/brachy_xcom_1.5MeV.muendat\n')
    egsinp.write('    muen for media = MUSCLE2_WW86, ADIPOSE2_WW86, SKIN2_WW86, GLAND2_WW86, CARTILAGE_WW86, HEART_BLOODFILLED_WW86, LUNG_BLOODFILLED_WW86, CORTICAL_BONE_WW86, YELLOW_MARROW_WW86, RED_MARROW_WW86, CALCIFICATION_ICRU46\n')
    egsinp.write('    dose scaling factor = %f\n' % DSF)
    egsinp.write(':stop scoring options:\n\n')

    egsinp.write('#' + '-' * 100 + '\n\n')
    egsinp.write('include file = /data/data062/sdeering/EGSnrc_Sept2016/egs_brachy/lib/transport/low_energy_default\n')

    egsinp.close()
    print('Created Egsinp File')

    if not os.path.exists('./Patient_Information/' + Name_String):
        os.makedirs('./Patient_Information/' + Name_String)
        print('Making Info Folder:' + './Patient_Information/' + Name_String)

    inputpath = './Patient_Information/%s/%s_Patient_%s_Info.txt' % (Name_String, Name_String, str(PatientNumber))
    info_file = open(inputpath, 'wb')
    info_file.write('This is the information for Patient %s %s File\n' % (str(PatientNumber),Name_String))
    info_file.write('File Created: %s \n\n' % str(datetime.date.today()))
    info_file.write('Boundaries and Cutoffs\n')
    info_file.write('------------------------------------------\n')
    info_file.write('Adipose/Gland Boundary Density       = %f\n' % AG_Cutoff_Den)
    info_file.write('Adipose/Gland Boundary HU            = %f\n' % AG_Cutoff_HU)
    info_file.write('Bone/Soft Tissue Boundary Density    = %f\n' % Bone_Cutoff_Den)
    info_file.write('Bone/Soft Tissue Boundary HU         = %f\n' % Bone_Cutoff_HU)
    info_file.write('Calcification/Gland Boundary Density = %f\n' % Calc_Cutoff_Den)
    info_file.write('Calcification/Gland Boundary HU      = %f\n\n' % Calc_Cutoff_HU)

    info_file.write('Metallic Artifact Reduction:\n')
    info_file.write('------------------------------------------\n')
    info_file.write('Replacement Value = ' + str(Replace) + '\n')
    info_file.write('Cutoff Value = ' + str(Cutoff) + '\n')
    info_file.write('Radius of Replacement Cylinder = ' + str(CylMaxRadius) + '\n')
    info_file.write('Number of High Density Voxels Replaced = ' + str(MAR_counter) + '\n\n')
    info_file.write('Number of Low Density Voxels Replaced = ' + str(low_artifact_counter) + '\n\n')

    info_file.write('Phantom Composition Information\n')
    info_file.write('------------------------------------------\n')
    info_file.write('Soft Tissue Voxels         = %d\n' % soft_tissue_count)
    info_file.write('Skin Voxels                = %d\n' % skin_count)
    info_file.write('Gland Voxels               = %d\n' % gland_count)
    info_file.write('Adipose Voxels             = %d\n' % adipose_count)
    info_file.write('Bone Voxels                = %d\n' % bone_count)
    info_file.write('Muscle Voxels              = %d\n' % muscle_count)
    info_file.write('Calcification Voxels       = %d\n' % calc_count)
    info_file.write('Air Voxels in Skin Contour = %d\n' % air_in_skin)
    info_file.write('Lead Covered Voxels        = 0\n')
    info_file.write('Outside of Skin Adipose Voxels Replaced: ' + str(adi_replace) + "\n\n")

    info_file.write('Structure List')
    info_file.write('------------------------------------------\n')
    for structures in structure_names:
        info_file.write(structures)
        info_file.write('  ')
    if roi_lung == -1:
        info_file.write('\n***********Warning: Could not find Lung Structure************')
    if roi_heart == -1:
        info_file.write('\n***********Warning: Could not find Heart Structure************')
    if roi_skin_margin == -1:
        info_file.write('\n***********Warning: Could not find Skin Margin Structure************')
    if roi_ptv == -1:
        info_file.write('\n***********Warning: Could not find PTV Structure************')
    if roi_ctv == -1:
        info_file.write('\n***********Warning: Could not find CTV Structure************')
    if roi_skin == -1:
        info_file.write('\n***********Warning: Could not find Skin Structure************')
    if roi_skin_surface == -1:
        info_file.write('\n***********Warning: Could not find Skin Surface Structure************')

    info_file.close()
    print('\n\n')

print('Script Finished: ' + str(sys.argv[0]))
print('END OF SCRIPT FOR PATIENTS:')
print(PatientNumbers)
print('Time Taken: ' + str(time.time() - start_time))

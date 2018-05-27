import sys
import dicom
import numpy as np
import os
from os import walk
import glob
import time
import copy
from scipy import spatial
from matplotlib.path import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pickle

def load_3ddose(PatientID, Label):  # Returns a tuple of the doses and errors, each in a 1D array.

    dose_filename = glob.glob('./3ddose/%s/*Pt_%s*' % (Label, PatientID))
    # print(dose_filename)
    c_dose_file = open(dose_filename[0], 'r+')

    c_data = []
    for line in c_dose_file:
        c_data.append(line)
    c_dose_file.close()
    print(len(c_data))
    dose_array = np.array(c_data[4].split(), dtype=float)
    error_array = np.array(c_data[5].split(), dtype=float)
    # print(len(dose_array))
    return dose_array,error_array

PatientNumbers = raw_input('Please enter 2-Digit Patient Number(s) (Delimited by Spaces) #: ')
Naming_String = raw_input('Name of 3DDose Files to Use = ')
Naming_String_TG43 = raw_input('Name of TG43 3DDose Files to Use = ')
start_time = time.time()

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
    #print(FilePath)
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
    ctfilenamelist = []
    for dirname, dirnames, filenames in walk(CTFilePath):
        ctfilenamelist.extend(filenames)
        break
    
    sort_dc = []
    for ct_file in ctfilenamelist:  #Removes .mim files from list
        if 'dcm' not in ct_file:
            ctfilenamelist.remove(ct_file)
            break                             #I think this break needs to go if there's more than just .mim files
       
        sort_dc.append((dicom.read_file(CTFilePath + ct_file).SliceLocation,CTFilePath + ct_file))
    
    sort_dc.sort() #Sorts the DCM Files from Inferior (-) to Superior (+)
    OrderedCTList = []
    for slice_location, DIFile in sort_dc:
        OrderedCTList.append(DIFile)
         
    #Now have ordered list of CT Dicom files for specified patient.
    #In the SliceListArray it goes [z][y][x]
    
    SliceListArray = []
    
    for current_ct in OrderedCTList:
        ct = dicom.read_file(current_ct)
        SliceListArray.append(ct.pixel_array)

    
#######################################################################
# Open Plan Dicom Files
#######################################################################    
    
    for root, dirnames, filenames in os.walk(PlanFilePath):
        for filename in filenames:
            if 'RTPLAN' in os.path.join(root,filename) and 'dcm' in os.path.join(root,filename):
                current_plan_file = os.path.join(root,filename)
    
    rp = dicom.read_file(current_plan_file)
    print('Finding Treatment Plan Information for Patient %s' % str(rp.PatientID))
                    
    # Extract the Seed Locations
    try:
        SeedLocations = []
        for seednumber,seedinfo in enumerate(rp.ApplicationSetupSequence):
            SeedLocation = []
    
            for i in range(0,3):
                SeedLocation.append(float(seedinfo.ChannelSequence[0].BrachyControlPointSequence[1].ControlPoint3DPosition[i]))
            
            # print(seednumber+1,SeedLocation)
            SeedLocations.append(SeedLocation) #Final List of Seed Locations
        
        print('Seed Locations Found')
    except:
        print('No Seed Locations Found for Patient # %d') % int(rp.PatientID)    
    
    # Get Essential Voxel and Slice Information from Reference Slice (First Slice)/
    
    Ref = dicom.read_file(OrderedCTList[0])
    Start_of_Grid = (float(Ref.ImagePositionPatient[0]), float(Ref.ImagePositionPatient[1]), float(Ref.SliceLocation))
    Size_of_Grid = (int(Ref.Columns), int(Ref.Rows), len(OrderedCTList))
    R_SOG = (Size_of_Grid[2], Size_of_Grid[1], Size_of_Grid[0])
    Voxel_Size = (float(Ref.PixelSpacing[0]),float(Ref.PixelSpacing[1]),2.00)
    Intercept = float(Ref.RescaleIntercept)
    Scaling = float(Ref.RescaleSlope)

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
            voxel_centers.append(
                ((xbounds[j] + 0.5 * Voxel_Size[0] / 10), (ybounds[i] + 0.5 * Voxel_Size[1] / 10)))
            # print(voxel_centers[len(voxel_centers)-513])

    ###################################################################
    ## Metallic Artifact Reduction
    ###################################################################

####################################################
## Write CT Images
###################################################################################

    MC_1D_Dose_Array_1 = load_3ddose(PatientNumber, Naming_String)[0]
    MC_1D_Dose_Array_2 = load_3ddose(PatientNumber, Naming_String_TG43)[0]
   #
    MC_3D_Dose_1 = np.reshape(MC_1D_Dose_Array_1, R_SOG)
    MC_3D_Dose_2 = np.reshape(MC_1D_Dose_Array_2, R_SOG)
   #  #RD_3D_Dose = Dose_Array_RD[0] * Dose_Scaling
   #
    Diff = MC_3D_Dose_1 - MC_3D_Dose_2
   #
   # # PerDiff1 = 100*np.divide(Diff1,MC_3D_Dose_1)
   #  #PerDiff2 = 100*np.divide(Diff2,MC_3D_Dose_2)

    x = np.asarray(xbounds)
    y = np.asarray(ybounds)

    for i in range(60,80):
        #sys.stdout.write('\r' + "Working on Slice " + str(i+1) + "/" + str(len(OrderedCTList)))
        #$sys.stdout.flush()    #Just for neatness in Terminal

        Slice_Seed_x = []
        Slice_Seed_y = []

        for Seeds in SeedLocations:
            sliceheight = int((Seeds[2] - 10 * zbounds[0]) // Voxel_Size[2])
            if sliceheight == i or sliceheight == (i - 1) or sliceheight == (i + 1):

                Slice_Seed_x.append(Seeds[0] / 10)
                Slice_Seed_y.append(Seeds[1] / 10)

        plt.pcolormesh(x, y, Diff[i], cmap='jet', vmin = -50, vmax = 50)
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        plt.colorbar()
        plt.text(-10, -20, str('Patient ' + PatientNumber + ' CT Slice ' + str(i)), color='white')
        # plt.show()
        plt.savefig('./COMP_Figures/CT_Dose_Diff_Pt_%s_Slice_%d.png' % (PatientNumber, i + 1), bbox_inches='tight', pad_inches=0)
        plt.close()

        #
        # plt.pcolormesh(x, y, MC_3D_Dose_2[i], vmin=0, vmax=250, cmap='jet')
        # plt.axis([x.min(), x.max(), y.min(), y.max()])
        # plt.colorbar()
        # plt.text(-20, -20, str(Naming_String_2), color='white')
        # # plt.show()
        # plt.savefig('./Dose_Plots/%s_vs_%s_Images/%s/%s_Dose_Slice_%s' % (Naming_String, Naming_String_2, PatientNumber, Naming_String_2, i + 1),bbox_inches='tight', pad_inches=0)
        # plt.close()
        #
        # plt.pcolormesh(x, y, Diff[i], vmin = -50, vmax = 50,  cmap='jet')
        # plt.axis([x.min(), x.max(), y.min(), y.max()])
        # plt.colorbar()
        # plt.text(-20, -20, str(i) + 'Difference', color='white')
        # #plt.show()
        # plt.savefig('./Dose_Plots/%s_vs_%s_Images/%s/Dose_Difference_%s' % (Naming_String, Naming_String_2, PatientNumber, i + 1),bbox_inches='tight', pad_inches=0)
        # plt.close()






    #print('....finished plotting')
    #print(CTFilePath)
    #print(PlanFilePath)
    #print(ContourFilePath)
      
#print('\nWrote Images to: ' + savefile)
print(PatientNumbers)
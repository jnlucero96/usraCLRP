import os
import time
import dicom
import glob
import numpy as np
import matplotlib.pyplot as plt
from os import walk
import pickle
import sys
import copy
from mpl_toolkits.mplot3d import Axes3D



def load_3ddose(PatientID, Label):  # Returns a tuple of the doses and errors, each in a 1D array.

    dose_filename = glob.glob('./3ddose/%s/*Pt_%s*' % (Label, PatientID))
    # print(dose_filename)
    c_dose_file = open(dose_filename[0], 'r+')

    c_data = []
    for line in c_dose_file:
        c_data.append(line)
    c_dose_file.close()
    #print(len(c_data))
    dose_array = np.array(c_data[4].split(), dtype=float)
    error_array = np.array(c_data[5].split(), dtype=float)
    # print(len(dose_array))
    return dose_array,error_array


def load_phant(patientID, label, p_width):

    phant_filename = glob.glob('./EGSphants/%s/*Pt_%s*' % (label, patientID))

    c_phant_file = open(phant_filename[0], 'r+')
    p_data = []
    comp_data = []
    for line in c_phant_file:
        p_data.append(line)
    c_phant_file.close()

    Number_of_Tissues = int(p_data[0])
    Tissue_Names = [] # !!! Tissue Names correspond exactly with numbers in comp_data. 0 = PlaceHolder.

    #print(p_data[0], p_data[5], p_data[22], p_data[180])

    for tissue_index in range(Number_of_Tissues):
        Tissue_Names.append(str.rstrip(p_data[tissue_index + 1], '\n'))

    print('Loading Phantom Information: ')

    for line in p_data:

        p_line = str.rstrip(line, '\n')
        # print(p_line)
        if len(p_line) == p_width:  # Tissue Assignment lines in EGSphants are one big string the size of the phantom width
            for digit in p_line:
                comp_data.append(int(digit))
    #print(len(comp_data))
    #print(comp_data[7])

    return Tissue_Names, np.array(comp_data)  # Tissue Names = List of Names,

def getKey3(item):
    return item[3]



#######################################################################################################################
# Accepting User Patient Input (Should have used Glob I think)
#######################################################################################################################

PatientNumbers = raw_input('Please enter 2-Digit Patient Number(s) (Delimited by Spaces) #: ')
conts_chosen = raw_input('Which Contour(s)(Delimited by Spaces): lung, heart, skin_margin, ptv1, ptv05, ctv, skin_skin, skin_surface, ribs, chest_wall, breast, body : ')
conts_chosen = conts_chosen.split(' ')

time.sleep(15000) #############TAKE OUT TO USEEEEEEEEE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    ###############################
    ## Get Patient CT Information
    ###############################

    ############################################################################
    # Opening and Sorting and Arraying the CT Files
    ############################################################################

    FirstFilePath = './MC PBSI Complete/Pt %s/' % PatientNumber
    FilePath = os.path.join(FirstFilePath, os.listdir(FirstFilePath)[0]) + '/'
    # print(FilePath)
    dir_sort = os.listdir(FilePath)
    for current_dir in dir_sort:
        if '_CT_' in current_dir:
            CTFilePath = FilePath + current_dir + '/'
        elif '_RTPLAN_' in current_dir:
            PlanFilePath = FilePath + current_dir + '/'
        elif '_RTst_' in current_dir:
            ContourFilePath = FilePath + current_dir + '/'

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

    #######################################################################
    # Open Plan Dicom Files
    #######################################################################

    for root, dirnames, filenames in os.walk(PlanFilePath):
        for filename in filenames:
            if 'RTPLAN' in os.path.join(root, filename) and 'dcm' in os.path.join(root, filename):
                current_plan_file = os.path.join(root, filename)

    rp = dicom.read_file(current_plan_file)
    print('Finding Treatment Plan Information for Patient %s' % str(rp.PatientID))

    # Extract the Seed Locations
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

    Ref = dicom.read_file(OrderedCTList[0])
    Ref2 = dicom.read_file(OrderedCTList[1])
    Start_of_Grid = (float(Ref.ImagePositionPatient[0]), float(Ref.ImagePositionPatient[1]), float(Ref.SliceLocation))
    Size_of_Grid = (int(Ref.Columns), int(Ref.Rows), len(OrderedCTList)) #(x,y,z)
    R_SOG = (Size_of_Grid[2],Size_of_Grid[1],Size_of_Grid[0]) # (z,y,x) Needed order for reshape I think
    Slice_Thickness = Ref2.SliceLocation - Ref.SliceLocation
    Voxel_Size = (float(Ref.PixelSpacing[0]), float(Ref.PixelSpacing[1]), Slice_Thickness)
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
        zbounds.append(round(zbounds[i] + Voxel_Size[2] / 10, 4))

    voxel_centers = []
    for i in range(len(ybounds) - 1):
        for j in range(len(xbounds) - 1):
            voxel_centers.append(((xbounds[j] + 0.5 * Voxel_Size[0] / 10), (ybounds[i] + 0.5 * Voxel_Size[1] / 10)))

    for working_contour in conts_chosen:

        if not os.path.exists(FilePath + '/Images/' + working_contour):
            os.makedirs(FilePath + '/Images/' + working_contour)
            print('Making Inputs Folder:' + FilePath + '/Images/' + working_contour)

        SliceListArray = []
        cont_array = []
        contmap = []

        for current_ct in OrderedCTList:
            #print('Starting a New Slice')
            ct = dicom.read_file(current_ct)
            SliceListArray.append(ct.pixel_array)

        cont_array = copy.deepcopy(SliceListArray)

        if os.path.isfile(FilePath + working_contour + '_contour.txt'):
            # print(FilePath + working_contour + '_contour.txt')
            contmap = pickle.load(open(FilePath + working_contour + '_contour.txt', 'rb'))
        else:
            print('Could not find ' + working_contour + ' contour')

        for z in range(R_SOG[0]):
            for y in range(R_SOG[1]):
                for x in range(R_SOG[2]):

                    if contmap[z][x + y * Size_of_Grid[0]]:
                        cont_array[z][y][x] = 6000

        print('Attempting to Write CT Images')
        if not os.path.exists(FilePath + 'Images/'):
            os.makedirs(FilePath + 'Images/')

        for i in range(len(SliceListArray)):

            Plot_This = False
            Slice_Seed_x = []
            Slice_Seed_y = []

            for Seeds in SeedLocations:
                sliceheight = int((Seeds[2] - 10 * zbounds[0]) // Voxel_Size[2])
                if sliceheight == i or sliceheight == (i - 1) or sliceheight == (i + 1):
                    Slice_Seed_x.append(Seeds[0] / 10)
                    Slice_Seed_y.append(Seeds[1] / 10)
                    Plot_This = True

            if Plot_This:

                x = np.asarray(xbounds)
                y = np.asarray(ybounds)
                plt.pcolormesh(x, y, cont_array[i], cmap='bone')
                plt.axis([x.min(), x.max(), y.min(), y.max()])
                plt.colorbar()
                plt.scatter(x=Slice_Seed_x, y=Slice_Seed_y, c='r', s=1, marker='x')
            #        plt.scatter(x=Skin_Margin_xs, y=Skin_Margin_ys, c='b', s=1, marker='-')
                plt.text(-20, -20, working_contour + ' contour', color='white')
                #plt.show()
                plt.savefig(FilePath + 'Images/%s/%s_contour_%03d.png' % (working_contour, working_contour, i + 1), bbox_inches='tight', pad_inches=0)
                plt.close()

    print('Finished with %s contour' % working_contour)

print('Script Finished: ' + str(sys.argv[0]))
print('END OF SCRIPT FOR PATIENTS:')
print(PatientNumbers)
print('Time Taken: ' + str(time.time() - start_time))

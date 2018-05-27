import os
import time
import dicom
import glob
import numpy as np
import matplotlib.pyplot as plt
from os import walk
import pickle


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



#######################################################################################################################
# Accepting User Patient Input (Should have used Glob I think)
#######################################################################################################################

PatientNumbers = raw_input('Please enter 2-Digit Patient Number(s) (Delimited by Spaces) #: ')
Naming_String = raw_input('Name of 3DDose Files to Use = ')
working_contours = raw_input('Contours to Use (Delimited by Spaces): ')
start_time = time.time()
write_trigger = False

if PatientNumbers == 'all':
    write_trigger = True
    pts = os.listdir('./MC PBSI Complete')
    PatientNumbers = []

    for i in pts:
        PatientNumbers.append(i.split('Pt ')[1])

else:
    PatientNumbers = PatientNumbers.split(' ')

working_contours = working_contours.split(' ')

d90_list = [[] for i in range(len(working_contours))]

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

    Ref = dicom.read_file(OrderedCTList[0])
    Size_of_Grid = (int(Ref.Columns), int(Ref.Rows), len(OrderedCTList)) #(x,y,z)
    R_SOG = (Size_of_Grid[2],Size_of_Grid[1],Size_of_Grid[0]) # (z,y,x) Needed order for reshape I think

    # Read in Monte Carlo Data from Both Phantoms Being Examined
    no_skip = True
    try:
        MC_1D_Dose_Array = load_3ddose(PatientNumber, Naming_String)[0]
        MC_Dose_Array = np.reshape(MC_1D_Dose_Array, R_SOG)
    except:
        no_skip = False

    if no_skip:

        for i,working_cont in enumerate(working_contours):

            cont_dose_array = []
            cont_dose_array_full = []

            if os.path.isfile(FilePath + working_cont + '_contour.txt'):
                current_contmap = pickle.load(open(FilePath + working_cont + '_contour.txt', 'rb'))

            else:
                print(FilePath + working_cont + '_contour.txt not found')
                print('Breaking out of ' + working_cont + ' loop for Pt #' + PatientNumber)
                if (working_cont == 'ptv1') or (working_cont == 'ptv05'):
                    print('Using CTV instead of PTV')
                    current_contmap = pickle.load(open(FilePath + 'ctv_contour.txt', 'rb'))
                else:
                    break

            for z in range(R_SOG[0]):
                for y in range(R_SOG[1]):
                    for x in range(R_SOG[2]):

                        if current_contmap[z][x + y * Size_of_Grid[0]]:

                            if MC_Dose_Array[z][y][x] > 499:
                                cont_dose_array.append(499.0)
                            else:
                                cont_dose_array.append(MC_Dose_Array[z][y][x])

                            cont_dose_array_full.append(MC_Dose_Array[z][y][x])

            ############################################################################################################
            ## How to do D90
            ############################################################################################################

            dose_range = (0,500)
            dose_range_for_value = (0,1000)
            weights = []

            bins = (dose_range[1] - dose_range[0]) * 10
            bins_for_value = (dose_range_for_value[1] - dose_range_for_value[0]) * 25
            cont_hist_full = plt.hist(cont_dose_array_full, bins=bins_for_value, range=dose_range_for_value, normed=1, histtype='step',cumulative=-1)
            max = 100
            for index,value in enumerate(cont_hist_full[0]):
                if abs(0.90000 - value) < max:
                    max = abs(0.90000 - value)
                    d90 = index

            Cont_Name = working_cont.upper()
            print('Patient # ' + str(PatientNumber))
            print(Naming_String + ' D90 for ' + Cont_Name + '  = ' + str(cont_hist_full[1][d90]) + ' Gy')
            plt.close()

            fig, ax = plt.subplots(figsize=(12, 6))
            cont_hist = ax.hist(cont_dose_array, bins=bins, range=dose_range, normed=1, histtype='step', cumulative=-1, label='D90 in ' + Cont_Name + ' - ' + Naming_String)
            ax.hlines(y=0.900, xmin=0, xmax=cont_hist_full[1][d90], linewidth=1, color='r')
            ax.text(300, 0.4, Cont_Name + "D90 = " + str(cont_hist_full[1][d90]) + ' Gy', color='b', size='large')
            ax.grid(True)
            ax.legend(loc='best')
            ax.set_title('Dose Metrics in ' + Cont_Name + ' Tissue: Patient ' + str(PatientNumber))
            ax.set_xlabel('Dose')
            ax.set_ylabel('Percentage of Voxels')
            #plt.show()
            fig.savefig('./Dose_Metrics/d90/' + Naming_String + '_' + working_cont + '_D90_Pt_' + str(PatientNumber))
            d90_list[i].append(cont_hist_full[1][d90])
    else:
        for i, working_cont in enumerate(working_contours):
            d90_list[i].append(0)

if write_trigger:
    for i,cont in enumerate(working_contours):
        pickle.dump(d90_list[i], open('./Dose_Metrics/metric_info/' + Naming_String + '_'+ cont + '_contour_ci_list.txt', 'wb'))


write_trigger = True
if write_trigger:
    read_file = open('./Dose_Metrics/metric_info/' + Naming_String + '_d90_read.txt', 'wb')
    read_file.write(' D90 Information for ' + Naming_String + '\n')

    for i,cont in enumerate(working_contours):
        read_file.write('\n\n-------- Contour Name: ' + cont + '---------\n\n')
        for b,PatientNumber in enumerate(PatientNumbers):
            read_file.write('Pt #%s D90 = %s \n' % (PatientNumber, str(d90_list[i][b])))







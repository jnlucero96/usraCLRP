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
from numpy import array, histogram
from matplotlib.pyplot import subplots, hist
from pickle import load, dump

# load self-modules
from py3ddose import DoseFile as read_dose

def plot_D90(
    patient_number, label, target_dir, d90_list, working_contour, 
    contour_dose_array, contour_dose_array_full
    ):

    """
    Description:
    Compute the dose to 90 percent of the patient. Plot the results.

    Inputs:
    :param :
    :type :

    Outputs:
    :param :
    :type :
    """

    weights = []

    dose_range = (0, 500)
    dose_range_for_value = (0, 1000)

    bins = 10 * (dose_range[1] - dose_range[0])
    bins_for_value = 25 * (dose_range_for_value[1] - dose_range_for_value[0])

    contour_hist_full, __, __ = hist(
        contour_dose_array_full, bins=bins_for_value, 
        range=dose_range_for_value, normed=True, histtype='step', 
        cumulative=-1
    )

    max_value = 100

    for index, value in enumerate(contour_hist_full[0]):
        if abs(0.90000 - value) < max_value:
            max_value = abs(0.9000 - value)
            d90_index = index
    
    contour_name = working_contour.upper()
    print "Patient # %s" % patient_number
    print label + " D90 for " + contour_name + " = " + \
    str(contour_hist_full[1][d90_index] + " Gy")

    fig, ax = subplots(1, 1, figsize=(10,10))
    contour_hist = ax.hist(
        contour_dose_array, bins=bins, range=dose_range, normed=True, 
        histtype='step', cumulative=-1
    )
    ax.hlines(
        y=0.900, xmin=0, xmax=contour_hist_full[1][d90_index],  linewidth=1, 
        color='r'
        )
    ax.grid(True)
    ax.legend(loc=0)
    ax.set_title(
        "Dose Metrics in " + contour_name + " Tissue: Patient " + 
        str(patient_number), fontsize=20
        )
    ax.set_xlabel("Dose (Gy)",fontsize=20)
    ax.set_ylabel("Percentage of Voxels", fontsize=20)

    fig.tight_layout()
    fig.savefig(
        target_dir + "/DOSE_METRICS/d90/" + contour_name + "_D90_Pt_" + 
        str(patient_number)
        )

    d90_list[i].append(contour_hist_full[1][d90_index])

    
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

    dir_path = target_dir + '/Pt_%s/' % patient_number
    patient_file_path = join(dir_path,listdir(dir_path)[0]) + '/'

    for current_dir in listdir(patient_file_path):
        if '_CT_' in current_dir:
            ct_file_path = patient_file_path + current_dir + '/'
        elif '_RTPLAN_' in current_dir:
            plan_file_path = patient_file_path + current_dir + '/'
        elif '_RTst_' in current_dir:
            contour_file_path = patient_file_path + current_dir + '/'
    
    print "Now importing and arranging CT data"
    print "=" * 50

    for __, __, filenames in walk(ct_file_path):  
        ct_filename_list.extend(filenames)
        break  # limit search to first directory
    
    for ct_file in ct_filename_list: 
        # for loop to remove .mim files
        if 'dcm' not in ct_file:
            ct_filename_list.remove(ct_file)
            #JNL: not sure why this break is here.
            break  
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
    # order: (x,y,z)
    grid_size = (int(ref.Columns), int(ref.Rows), len(ordered_ct_list)) 
    # order: (z,y,x) needed order for reshape
    r_sog = (grid_size[2], grid_size[1], grid_size[0]) 

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

            contmap_file = patient_file_path + working_contour + '_contour.txt'

            if isfile(contmap_file):
                current_contmap = load(
                    open(contmap_file,'rb')
                )
            else:
                print contmap_file + " not found"
                print "Breaking out of " + working_contour + \
                " loop for Pt # " + patient_number
                if working_contour == 'ptv1' or working_contour == 'ptv05':
                    print "Using CTV instead of PTV"
                    current_contmap = load(
                        open(contmap_file, 'rb')
                    )
                else:
                    break
            
            for z_index in xrange(r_sog[0]):
                for y_index in xrange(r_sog[1]):
                    for x_index in xrange(r_sog[2]):
                        if current_contmap[z_index][
                            x_index + y_index * grid_size[0]
                            ]:

                            #JNL: not sure what the 499 is for in the check
                            if mc_dose_array[z_index][y_index][x_index] > 499: 
                                cont_dose_array.append(499.0)
                            else:
                                cont_dose_array.append(
                                    mc_dose_array[z_index][y_index][x_index]
                                )
                            
                            cont_dose_array_full.append(
                                mc_dose_array[z_index][y_index][x_index]
                            )

            plot_D90() #JNL-TODO: fill in the neccessary arguments
        
        else:
            for contour_index, working_contour in enumerate(working_contour):
                d90_list[i].append(0.0)

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
    
        get_CT() #JNL-TODO: Finish inputting the arguments into this function

    if write_trigger:
        for index, contour in enumerate(working_contours):
            with open(
                target_dir + '/Dose_Metrics/metric_info/' + naming_string + 
                '_' + contour + '_contour_ci_list.txt', 'wb'
                ) as d90_file:
            
                dump(d90_list[index], d90_file)
    
    with open(
        target_dir + '/Dose_Metrics/metric_info' + naming_string 
        + '_d90_read.txt', 'wb'
        ) as info_file:

        info_file.write(" D90 Information for " + naming_string + "\n")

        for index, contour in enumerate(working_contours):
            info_file.write(
                "\n\n ---------- Contour Name: " + contour 
                + "---------\n\n" 
                )

            for index2, patient_number in enumerate(patient_numbers):
                info_file.write(
                    "Pt #%s D90 = %s \n" % (
                        patient_number, d90_list[index][index2]
                        )
                )

    



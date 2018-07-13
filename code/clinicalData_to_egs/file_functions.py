#!/usr/bin/env
# author: Joseph Lucero
# created on: 25 May 2018 22:13:43
# purpose: define functions that handles files

# import external functions
from os import getcwd, makedirs
from os.path import exists, isdir
from glob import glob
from numpy import array, empty, interp

from get_funcs import get_target_dir, get_CT_calibration, get_media

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
            
def create_EGSPhant(
    slice_array_lst, ref_intercept, SIZE_OF_GRID, bounds, 
    path_to_calibration=None, path_to_egsphants=None,
    egsphant_name='my_egsphant'
    ):

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

    CT_HU, CT_rho = get_CT_calibration(path_to_calibration)

    [
        media_names, media_HU_ranges, media_density, PEGS4_names, HU_min, HU_max
        ] = get_media()

    xbounds, ybounds, zbounds = bounds

    if not path_to_egsphants:
        while True:
            path_to_egsphants = raw_input(
                "Please input the full directory path to where you would like the egsphants to go:>> "
            )
            if not isdir(path_to_egsphants):
                print \
                "This input is not a valid directory." 
                answer = raw_input(
                    "Would you like me to make that folder? [y/n]:>> "
                    )
                if answer.lower() in ('y','yes'):
                    print "Making egsphant directory as requested."
                    makedirs(path_to_egsphants)
            else:
                break

    egsphant_file_path = input_file_path = path_to_egsphants + \
        '/{0}.egsphant'.format(egsphant_name)

    with open(egsphant_file_path,'w+') as egsphant:
        egsphant.write(
            '{0}'.format(len(media_names)) + '\n'
        )
        for media in PEGS4_names:
            egsphant.write(
                media + '\n'
            )

        egsphant.write('  0.25')
        for __ in xrange(len(media_names) - 1):
            egsphant.write(' '*7 + '0.25')
        egsphant.write('\n')

        for dim_size in SIZE_OF_GRID:
            egsphant.write(' '*2 + '{0}'.format(dim_size))

        newline_counter = 0 
        print "Writing voxel boundaries..."
        for x_bound in xbounds:
            if newline_counter % 3 == 0:
                egsphant.write('\n' + '\t')

            egsphant.write("{0:.5f}".format(x_bound))

            if newline_counter % 3 !=2:
                egsphant.write(' ' * 10)
            elif newline_counter % 3 == 2:
                egsphant.write(' ' * 7)
            
            newline_counter += 1

        for y_bound in ybounds:
            if newline_counter % 3 == 0:
                egsphant.write('\n' + '\t')

            egsphant.write("{0:.5f}".format(y_bound))

            if newline_counter % 3 != 2:
                egsphant.write(' ' * 10)
            elif newline_counter % 3 == 2:
                egsphant.write(' ' * 7)

            newline_counter += 1

        for z_bounds in zbounds:
            if newline_counter % 3 == 0:
                egsphant.write('\n' + '\t')

            egsphant.write("{0:.5f}".format(z_bounds))

            if newline_counter % 3 != 2:
                egsphant.write(' ' * 10)
            elif newline_counter % 3 == 2:
                egsphant.write(' ' * 7)

            newline_counter += 1
        
        egsphant.write('\n')

        media_dict = {media:0 for media in media_names.upper()}

        print "Writing tisue assignments..."

        for k in xrange(SIZE_OF_GRID[2]):
            for j in xrange(SIZE_OF_GRID[1]):
                for i in xrange(SIZE_OF_GRID[0]):

                    ct_value = slice_array_lst[k][j][i] + ref_intercept 

                    for (
                        index, (media_HU_low_bound, media_HU_up_bound)
                        ) in enumerate(media_HU_ranges, 1):
                        if media_HU_low_bound <= ct_value <= media_HU_up_bound:
                            egsphant.write('{0}'.format(index)) 
                egsphant.write('\n')
            egsphant.write('\n')

        print "Finished with tissue assignments"

        print "Writing voxel density assignments..."

        newline_counter = 0  # reset newline counter

        for z in xrange(SIZE_OF_GRID[2]):
            print "\r" + str(z) + '/' + str(SIZE_OF_GRID[2])
            for y in xrange(SIZE_OF_GRID[1]):
                for x in xrange(SIZE_OF_GRID[0]):
                    ct_value = slice_array_lst[z][y][x] + ref_intercept
                    current_density = interp(ct_value, CT_HU, CT_rho)
                    if newline_counter % 3 == 0:
                        egsphant.write('\n' + '  ')
                    egsphant.write('{0:.15f}'.format(current_density))

                    if newline_counter % 3 != 2:
                        egsphant.write(' ' * 9)
                    elif newline_counter % 3 == 2:
                        egsphant.write(' ' * 7)
                    newline_counter += 1

        print "Finished writing voxel density assignments..."

def create_egsinp():

    """
    Description:

    Inputs:

    Outputs:

    """

    


        


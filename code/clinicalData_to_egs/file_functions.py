#!/usr/bin/env python3
# author: Joseph Lucero
# created on: 25 May 2018 22:13:43
# purpose: define functions that handles files

# import external functions
from datetime import date
from os import getcwd, makedirs
from os.path import exists, isdir, expanduser
from glob import glob
from numpy import array, empty, interp

from get_funcs import (
    get_target_dir, get_CT_calibration, get_media, get_egs_brachy_settings
    )

def load_phantom(patientID, label, p_width):
    """\
    Description: 
    Function that takes in patient ID and label to return the 
    doses and errors associated with that patient in 1D arrays.

    *** currently not incorporated into main code ***\
    """

    p_data = []
    comp_data = []
    tissue_names = empty(p_width, dtype=int)
    target_dir = get_target_dir()

    while True:
        phantom_filename = input("Please input path to target file: ")
        if not exists(target_dir):
            print("That was not a valid path.")
            continue
        else:
            break
    
    with open(phantom_filename,'r+') as c_phantom_file:
        for line in c_phantom_file:
            p_data.append(line)
    
    num_of_tissues = int(p_data)
    for tissue_index in range(num_of_tissues):
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
    slice_array_lst, contour_map, ref_intercept, SIZE_OF_GRID, bounds, 
    path_to_calibration=None, path_to_media=None, path_to_egsphants=None,
    egsphant_name='my_egsphant'
    ):

    """\
    Description: 
    Script to create EGSPhant files from CT files given a specific 
    CT calibration curve.

    Inputs:
    :param slice_array_lst: Pixel intensity data in every slice of the CT image
    :type slice_array_lst: list
    :param contour_map: List of all the contours in the image and the shape \
    of the contours
    :type contour_map: list
    :param ref_intercept: Rescaling constant for the HU readings
    :type ref_intercept: int
    :param SIZE_OF_GRID: Number of voxels in each given dimension
    :type SIZE_OF_GRID: tuple
    :param bounds: Coordinates of the bounding planes in each of the 3 \
    dimensions
    :type bounds: list
    :param path_to_calibration: Full file path to where the CT calibration \
    file is stored
    :type path_to_calibration: str
    :param path_to_media: Full file path to where the tissue assignment file \
    is stored
    :type path_to_media: str
    :param path_to_egsphants: Full directory path to where the egsphants are \
    to be stored after they are created
    :type path_to_egsphants: str
    :param egsphant_name: The BASENAME of the egsphant so the output file has \
    the name BASENAME.egsphant
    :type egsphant_name: str

    Outputs:
    :param egsphant: egsphant file created in the specified directory
    :type egsphant: egsphant file\
    """

    print("\n\n" + "="*50)
    print("Attempting to create egsphant files...")
    print("="*50)

    CT_HU, CT_rho = get_CT_calibration(path_to_calibration)

    [
        media_names, media_HU_ranges, media_HU, media_density, PEGS4_names
        ] = get_media(path_to_media)

    xbounds, ybounds, zbounds = bounds

    if not path_to_egsphants:
        while True:
            path_to_egsphants = expanduser(input(
                "Please input the full directory path to where you would like the egsphants to go:>> "
            ))
            if not isdir(path_to_egsphants):
                print("This input is not a valid directory.")
                answer = input("Would you like me to make that folder? [y/n]:>> ")
                if answer.lower() in ('y','yes'):
                    print("Making egsphant directory as requested.")
                    makedirs(expanduser(path_to_egsphants))
                    break
            else:
                break

    egsphant_file_path = path_to_egsphants + \
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
        for __ in range(len(media_names) - 1):
            egsphant.write(' '*7 + '0.25')
        egsphant.write('\n')

        for dim_size in SIZE_OF_GRID:
            egsphant.write(' '*2 + '{0}'.format(dim_size))

        newline_counter = 0 
        print("Writing voxel boundaries...")
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

        for z_bound in zbounds:
            if newline_counter % 3 == 0:
                egsphant.write('\n' + '\t')

            egsphant.write("{0:.5f}".format(z_bound))

            if newline_counter % 3 != 2:
                egsphant.write(' ' * 10)
            elif newline_counter % 3 == 2:
                egsphant.write(' ' * 7)

            newline_counter += 1
        
        egsphant.write('\n')
        egsphant.flush()

        media_dict = {media.upper(): 0 for media in media_names}

        print("Writing tisue assignments...")

        for k in range(SIZE_OF_GRID[2]):
            print("\r" + str(k) + '/' + str(SIZE_OF_GRID[2]))
            for j in range(SIZE_OF_GRID[1]):
                for i in range(SIZE_OF_GRID[0]):

                    ct_value = slice_array_lst[k][j][i] + ref_intercept 

                    # This for-loop is not robust if you have overlapping 
                    # boundaries for HU min and HU max of different media...
                    # Ask Rowan about how to improve it. Maybe add contour 
                    # information to break degeneracies?
                    for (
                        index, (media_HU_low_bound, media_HU_up_bound)
                        ) in enumerate(media_HU_ranges, 1):
                        if media_HU_low_bound <= ct_value <= media_HU_up_bound:
                            egsphant.write('{0}'.format(index)) 
                            media_dict[media_names[index-1].upper()] += 1
                            break;
                    else:
                        print("No media found for these set of indexes:", tuple(i,j,k))
                egsphant.write('\n')
            egsphant.write('\n')
        
        egsphant.flush()
        print("Finished with tissue assignments")
        print("Media Count:")
        print(media_dict)

        print()
        print("Writing voxel density assignments...")

        newline_counter = 0  # reset newline counter

        for z in range(SIZE_OF_GRID[2]):
            print("\r" + str(z) + '/' + str(SIZE_OF_GRID[2]))
            for y in range(SIZE_OF_GRID[1]):
                for x in range(SIZE_OF_GRID[0]):
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
                egsphant.write('\n')
                newline_counter = 0
            egsphant.write('\n')
        egsphant.flush()
        print("Finished writing voxel density assignments...")

    print("Finished writing egsphant...")

def create_egsinp(egsinp_path, script_name):

    """\
    Description: Given user defaults construct egsinp file that can be used to
    run egs_brachy simulations 

    *** currently not incorporated into main code ***\
    """

    print("\n Attempting to create input file")
    print('-'*50)

    format_dict, inp_settings = get_egs_brachy_settings(script_name)

    with open(egsinp_path, 'w') as egsinp:
        egsinp.write(
            """\
            ###################################################################
            ######
            ###### This is an automatically generated egsinp file created on {format_dict[date]}
            ###### Information for this file can be found at: {format_dict[info_dir]}
            ###### File created by script: {format_dict[file_source]}
            ######
            ###################################################################

            #------------------------------------------------------------------
            :start run control:
                ncase = {inp_settings[ncase]}
                nbatch = {inp_settings[nbatch]}
                nchunk = {inp_settings[nchunck]}
                calculation = {inp_settings[calculation]}
                geometry error limit = {inp_settings[geom_err_lim]}
            :stop run control:

            #------------------------------------------------------------------
            :start run mode:
                run mode = {inp_settings[run_mode]}
            :stop run mode: 

            #------------------------------------------------------------------
            :start media definition:
                AE = {inp_settings[AE]}
                UE = {inp_settings[UE]}
                AP = {inp_settings[AP]}
                UP = {inp_settings[UP]}

                material data file = {inp_settings[material_dat_file_path]}
            :stop media definition:

            #------------------------------------------------------------------
            :start geometry definition:
                
                :start geometry:
                    name = phantom
                    library = egs_glib
                    type = egsphant
                    egsphant file = {inp_settings[egsphant_file_path]}
                    density file = {inp_settings[density_file_path]}
                :stop geometry:

                :start geometry:
                    name = seed
                    library = egs_glib
                    include file = {inp_settings[seed_file_path]}
                :stop geometry:

                :start geometry:
                    name = phantom_w_seeds
                    library = egs_autoenvelope
                    :start inscribed geometry:
                        inscribed geometry name = seed
                        
                        :start transformations:
                            include file = {inp_settings[transform_file_path]}
                        :stop transformations:
                        
                        :start volume correction:
                            correction type = correct
                            density of random points (cm^-3) = {inp_settings[volcor_density]}
                            include file = {inp_settings[seed_boundary_file_path]}
                        :stop volume correction
                    :stop inscribed geometry:
                :stop geometry:
            :stop geometry definition:

            #------------------------------------------------------------------
            :start volume correction:
                :start source volume correction:
                    correction type = {inp_settings[volcor_type]}
                    density of random points (cm^-3) = 1E8
                    include file = {inp_settings[seed_boundary_file_path]}
                :stop source volume correction:
            :stop volume correction:

            #------------------------------------------------------------------
            :start source definition:
                :start source:
                    name = source_name
                    library = egs_isotropic_source
                    charge = 0
                    include file = {inp_settings[seed_shape_file_path]}

                    :start spectrum:
                        type = tabulated spectrum
                        spectrum file = {inp_settings[spectrum_file_path]}
                    :stop spectrum:
                :stop source:

                :start transformations:
                    include file = {inp_settings[transform_file_path]}
                :stop transformations:

                simulation source = source_name
            :stop source defintion:

            #------------------------------------------------------------------
            :start scoring options:
                score tracklength dose = {inp_settings[score_tracklength]}
                score energy deposition = {inp_settings[score_deposition]}
                score scatter dose = {inp_settings[score_scatter]}
                muen file = {inp_settings[muen_file_path]}
                muen for media = {inp_settings[muen_for_media]}
                
                dose file format = gzip

                output egsphant files = {inp_settings[output_egsphant]}
                egsphant file format = gzip

                output voxel info files = {inp_settings[output_voxel_info]}
                voxel info file format = gzip

                output volume correction file for phantoms = {inp_settings[output_volcor_file]}
                volume correction file format = gzip

                record initial position = {inp_settings[record_init_pos]}

                dose scaling factor = {inp_settings[dose_scaling_factor]}
            :stop scoring options:

            #------------------------------------------------------------------
            # Transport parameters
            include file = {inp_settings[tranport_param_file_path]} \
            """.format(format_dict=format_dict, inp_settings=inp_settings)
        )


    


        


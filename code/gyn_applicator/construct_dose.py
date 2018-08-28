#!/usr/bin/env python3

from gzip import open as gOpen
from glob import glob
from os.path import expanduser, isfile, isdir
from os import getcwd
import numpy

# define dose file object 
# modified from github/christopherpoole/3DDose
class DoseFile(object):

    """Description: Loads 3ddose file into numpy arrays"""

    def __init__(self, file_name, load_uncertainty=False):
        self._load_3ddose(file_name, load_uncertainty)
        self.shape = self.dose.shape
        self.size = self.dose.size
    
    def _load_3ddose(self, file_name, load_uncertainty=False):

        if file_name[-3:] == '.gz':  # open gzipped file
            with gOpen(file_name, 'rb') as data_file:
                data = data_file.read().split('\n')
        else:
            with open(file_name, 'rb') as data_file: # open regular 3ddose file
                data = data_file.read().split('\n')
        
        cur_line = 0
        
        x, y, z = map(int, data[cur_line].split())
        self.shape = (z,y,x)
        # print self.shape
        self.size = numpy.multiply.reduce(self.shape)
        
        cur_line += 1
        
        positions = []
        for i in range(0,3):
            bounds = []
            while len(bounds) < [x, y, z][i]:
                line_positions = map(float, data[cur_line].split())
                bounds += line_positions
                cur_line += 1
            positions.append(bounds)
        
        # recall that dimensions are read in order x, y, z
        positions = [positions[0], positions[1], positions[2]]
        
        self.positions = positions
        
        dose = []
        while len(dose) < self.size:
            line_data = map(float, data[cur_line].split())
            dose += line_data
            cur_line += 1
        
        # have dose array be in order (x, y, z)
        self.dose = self.dose.reshape((self.shape)).transpose((2,1,0))
        
        if load_uncertainty:
            uncertainty = []
            while len(uncertainty) < self.size:
                line_data = map(float, data[cur_line].split())
                uncertainty += line_data
                cur_line += 1
            self.uncertainty = numpy.array(uncertainty)  
            assert len(self.uncertainty) == self.size, \
            "len of uncertainty = {} (expected {})".format(len(self.uncertainty), self.size)
            # have uncertainty array be in order (x, y, z)
            self.uncertainty = self.uncertainty.reshape((self.shape)).transpose((2,1,0))
            assert self.uncertainty.size == self.size, \
            "uncertainty array size does not match that specified."

        # reset the shape of array to have the proper (x,y,z) ordering
        self.shape = (x, y, z)

def checks(weights, len_weight, len_dose_files):

    """\
    Description: Does a few checks that happen before the weighted superposition \
    of the reference dose files is created. Specifically checks for whether \
    enough weights have been specified for the reference dose files. 

    Inputs:
    :param weights: series of weights that are used during the superposition
    :type weights: numpy.array or list
    :param len_weight: measure of how many weights have been specified
    :type len_weight: numpy.array or list
    :param len_dose_files: measure of how many reference files there are 
    :type len_dose_files: int

    Outputs:
    :param weights: series of weights, where number of weights matches the \
    number of reference files
    :type weights: numpy.ndarray or list
    :param type_of_weights: tells whether weights array is a numpy array or \
    Python list
    :type type_of_weights: type\
    """

    type_of_weights = type(weights)

    if len_weight != len_dose_files:

        print("WARNING: Number of weights specified not the same as the number \
        of reference dose files specified.") 

        difference = len_dose_files - len_weight

        print("Should we set the unspecified weights to unity [y/n]?")
        while True:
            
            answer = input(">> ")

            if answer.lower() in ('y', 'yes'):
                
                if type_of_weights is numpy.ndarray:
                    array_to_append = numpy.ones(len_dose_files - len_weight)
                    weights = numpy.append(weights, array_to_append)
                elif type_of_weights is list:
                    
                    counter = 0
                    while counter < difference:
                        weights.append(1.0)
                        counter += 1
                else:
                    print(
                        "ALERT: Weights should either be array or list. Type \
                        not understood. Exiting program now..."
                        ); exit(1)
                break
            elif answer.lower() in ('n', 'no'):
                print("Please provide the remaining weights in space-delimited format")
                remaining_list = list(map(float, input(">> ").split()))
                
                while len(remaining_list) != abs(difference):
                    print(
                        "You have not specified the correct number of \
                        remaining weights. The number of remaining weights \
                        that you need to specify is {0}. Please input the \
                        weights again.".format(difference)
                        )
                    remaining_list = list(map(float, input(">> ").split()))

                if type_of_weights is numpy.ndarray:
                    array_to_append = numpy.array(remaining_list)
                    weights = numpy.append(weights, array_to_append)
                elif type_of_weights is list:
                    for remaining_weight in remaining_list:
                        weights.append(remaining_weight)
                else:
                    print(
                        "ALERT: Weights should either be array or list. Type \
                        not understood. Exiting program now..."
                        ); exit(1)
                break
            else: 
                print("Input not understood. Please input answer again [y/n].")
                continue
                
    else:
        pass

    return weights, type_of_weights

def construct_dose(weights_input, ref_dose_files_input, get_uncertainty):

    """\
    Description: Construct a dose profile from the reference dose files using \
    linear superposition where the given weights are specified by the user. 

    Inputs:
    :param weights: weights with which the reference dose files are to be added.
    :type weights: array or list of floats
    :param dose_files: reference dose files that are to be added
    :type dose_files: list of strings
    :param get_uncertainty: output uncertainty or not
    :type dose_files: bool

    Output: 
    :param dose_array: the dose array created by the weighted superposition of \
    the reference dose array
    :type dose_array: numpy.ndarray
    :param resultant_uncertainty: the absolute uncertainty of the dose array
    :type resultant_uncertainty: numpy.ndarray\
    """

    weights_unscaled, type_of_weights = checks(
        weights_input, len(weights_input), len(ref_dose_files_input)
        )

    max_weight = max(weights_input)

    if type_of_weights is numpy.ndarray:
        weights = weights_unscaled / max_weight
    elif type_of_weights is list:
        weights = [weight / max_weight for weight in weights_unscaled]

    for index, ref_file in enumerate(ref_dose_files_input):
        data = DoseFile(ref_file, load_uncertainty=get_uncertainty)
        if index == 0:
            resultant_dose = resultant_uncertainty = numpy.zeros(data.shape)
        
        resultant_dose += weights[index] * data.dose
        if get_uncertainty:
            resultant_uncertainty += (
                weights[index] * data.uncertainty * data.dose
                ) ** 2

    resultant_uncertainty = numpy.sqrt(resultant_uncertainty)

    x_pos, y_pos, z_pos, = numpy.array(data.positions)
    x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2
    y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2
    z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2

    position_tuple = (x_pos_mid, y_pos_mid, z_pos_mid)

    return resultant_dose, resultant_uncertainty, position_tuple
        
def main():

    """\
    Description: Run superposition calculation.

    Outputs:
    :param gfile: Output dat file containing the resultant dose from doing \
    weighted superposition of the reference dose files. Has 5 columns containing, \
    in order, the x-coordinate, the y-coordinate, the z-coordinate, the dose at coordinate, \
    and error of dose at coordinate.
    :type gfile: gzipped data file\
    """

    cwd = getcwd()

    print("Please input the size of the Nucletron applicator you wish to work with.")
    applicator = input(">> ")
    if applicator not in ('25', '30', '35', '40'):
        while applicator not in ('25', '30', '35', '40'):
            print("That is not a valid applicator size. Please input in another size.")
            applicator = input(">> ")

    target_dir = cwd + '/reference/{0}mm_applicator'.format(applicator)
    
    if not isdir(target_dir):
        print("Cannot find the reference directory containing the reference dose files.")
        print("Please ensure that the reference directory is in the same directory as this script.")
        exit(1)
    
    ref_dose_files = [ref_file for ref_file in glob(target_dir + '/*.3ddose.gz')]
    print("Do you want to use the default weights [y/n]?")
    while True:
        answer = input(">> ")
        if answer.lower() in ('y', 'yes'):
            print("Setting the default weights.")
            weights_input = numpy.array([
                0.32032854, 0.5687885, 1.00, 0.91991786, 0.49897331,
                0.3100616, 0.32238193, 0.3798768, 0.59958932, 0.82956879, 1.00
                ]) # standard weights from Ali & Cygler
            break
        elif answer.lower() in ('n', 'no'):
            print(
                "Please input weights for the reference dose files separated by \
                whitespace. You need to input "
                + "in {0} weights".format(len(ref_dose_files))
            )
            user_input = input(">> ")
            if isfile(expanduser(user_input)):
                weights_input = numpy.loadtxt(expanduser(user_input))
            else:
                weights_input = list(map(float, user_input.split()))
            break
        else:
            print "Input not understood."

    # construct the dose from weighted superposition of the reference dose files
    dose, uncertainty, positions = construct_dose(
        weights_input, ref_dose_files, get_uncertainty=True
        )

    # write into file
    with gOpen('resultant_file.dat', 'w') as gfile:
        for x_index, x in positions[0]:
            for y_index, y in positions[1]:
                for z_index, z in positions[2]:
                    gfile.write(
                        '{0}\t{1}\t{2}\t{3:.15f}\t{4:.15f}\n'.format(
                            x, y, z, 
                            dose[x_index, y_index, z_index], 
                            uncertainty[x_index, y_index, z_index]
                            )
                        )   
                
    return 0

if  __name__ == "__main__":
    main()


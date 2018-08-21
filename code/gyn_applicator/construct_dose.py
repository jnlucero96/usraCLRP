#!/usr/bin/env python3

from gzip import open as gOpen
from glob import glob
from os.path import expanduser
import numpy

# define dose file object 
class DoseFile(object):

    # copied from github/christopherpoole/3DDose

    def __init__(self, file_name, load_uncertainty=False):
        """
        Attempts to detect the dose file etension automatically. If an unknown
        extension is detected, loads a .3ddose file by default.
        """

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
        self.spacing = [numpy.diff(p) for p in self.positions]
        self.resolution = [s[0] for s in self.spacing if s.all()]
        self.origin = numpy.add( 
            [p[0] for p in positions], numpy.array(self.resolution)/2.
            )
        
        assert len(self.resolution) == 3, \
        "Non-linear resolution in either x, y or z."
        
        dose = []
        while len(dose) < self.size:
            line_data = map(float, data[cur_line].split())
            dose += line_data
            cur_line += 1
        
        self.dose = numpy.array(dose)
        assert len(self.dose) == self.size, \
        "len of dose = {} (expect {})".format(len(self.dose), self.size)
        # have dose array be in order (x, y, z)
        self.dose = self.dose.reshape((self.shape)).transpose((2,1,0))
        assert self.dose.size == self.size, \
        "Dose array size does not match that specified."
        
        if load_uncertainty:
            uncertainty = []
            while len(uncertainty) < self.size:
                line_data = map(float, data[cur_line].split())
                uncertainty += line_data
                cur_line += 1
            self.uncertainty = numpy.array(uncertainty)  
            assert len(self.uncertainty) == self.size, \
            "len of uncertainty = {} (expected {})".format(
                len(self.uncertainty), self.size
                )
            # have uncertainty array be in order (x, y, z)
            self.uncertainty = self.uncertainty.reshape(
                (self.shape)
                ).transpose((2,1,0))
            assert self.uncertainty.size == self.size, \
            "uncertainty array size does not match that specified."

        # reset the shape of array to have the proper (x,y,z) ordering
        self.shape = (x, y, z)

    def max(self):
        return self.dose.max()
        
    def min(self):
        return self.dose.min()

    @property
    def x_extent(self):
        return self.positions[0][0], self.positions[0][-1]

    @property
    def y_extent(self):
        return self.positions[1][0], self.positions[1][-1]
    
    @property
    def z_extent(self):
        return self.positions[2][0], self.positions[2][-1]

def checks(weights, len_weight, ref_dose_files, len_dose_files):

    type_of_weights = type(weights)

    if len_weight != len_dose_files:

        difference = len_dose_files - len_weight

        if difference > 0:
            print("WARNING: Number of weights specified not the same as the number \
            of reference dose files specified.") 

            while True:
                print("Should we set the unspecified weights to unity [y/n]?")
                answer = input(">> ")

                if answer.lower() in ('y', 'yes'):
                    
                    if type_of_weights is numpy.ndarray:
                        array_to_append = numpy.ones(len_dose_files - len_weight)
                        numpy.append(weights, array_to_append)
                    elif type_of_weights is list:
                        difference = len_dose_files - len_weight
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
                        numpy.append(weights, array_to_append)
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
                    
        elif difference < 0:
            print("WARNING: Number of reference dose files specified not the same \
            as the number weights specified.")

            while True: 
                print("Do you want to ignore the extra specified weights")
                answer = input(">> ")
                if answer in ('y', 'yes'):
                    print("Ignoring the extra specified weights.")
                    weights = weights[:len_dose_files]
                elif answer in ('n', 'no'): 
                    print("Please input the correct weights.")
                    correct_weights = list(map(float, input(">> ").split()))
                    
                    while len(correct_weights) != len_dose_files:
                        print(
                            "You have not specified the correct number of \
                            remaining weights. The number of weights \
                            that you need to specify is {0}. Please input the \
                            weights again.".format(len_dose_files)
                            )
                        correct_weights = list(map(float, input(">> ").split()))

                    if type_of_weights is numpy.ndarray:
                        weights = numpy.array(correct_weights)
                    elif type_of_weights is list:
                        weights = correct_weights
                    else:
                        print(
                            "ALERT: Weights should either be array or list. Type \
                            not understood. Exiting program now..."
                            ); exit(1)
                    break
                else: 
                    print(
                        "Input not understood. Please input answer again \
                        [y/n]."
                        )
                    continue
    else:
        pass

    return weights, ref_dose_files, type_of_weights

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
    :param dose_array:
    :type dose_array:
    :param resultant_uncertainty:
    :type resultant_uncertainty:
    """

    weights_unscaled, ref_dose_files, type_of_weights = checks(
        weights_input, len(weights_input), ref_dose_files_input, 
        len(ref_dose_files_input)
        )

    max_weight = max(weights)

    if type_of_weights is numpy.ndarray:
        weights = weights_unscaled / max_weight
    elif type_of_weights is list:
        weights = [weight / max_weight for weight in weights_unscaled]

    for index, ref_file in enumerate(ref_dose_files):
        data = DoseFile(ref_file, load_uncertainty=get_uncertainty)
        if index = 0:
            resultant_dose = resultant_uncertainty = zeros(data.shape)
        
        resultant_dose += weights[index] * data.dose
        if get_uncertainty:
            resultant_uncertainty += (weights[index] * data.uncertainty) ** 2

    resultant_uncertainty = numpy.sqrt(resultant_uncertainty)

    return resultant_dose, resultant_uncertainty
        
def main():

    print(
        "Please input the file path to the directory containing \
        reference dose files"
        )
    target_dir = expanduser(input(">> "))
    ref_dose_files = [ref_file for ref_file in glob(target_dir + '/*.3ddose')]
    print(
        "Please input weights for the reference dose files separated by \
        whitespace. You need to input in {0} weights".format(len(ref_dose_files))
        )
    weights_input = list(map(float, input(">> ").split()))

    dose, uncertainty = construct_dose(
        weights_input, ref_dose_files, get_uncertainty=
        )

    if save_dose:
        with gOpen('resultant_file.dat', 'w') as gfile:
            for dose_value, uncertainty_val in zip(dose, uncertainty):
                gfile.write(
                    '{0:.15f}\t {1:.15f}\n'.format(dose_value, uncertainty_val)
                    )

    

    return 0


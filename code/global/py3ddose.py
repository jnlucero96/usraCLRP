#!/usr/bin/env
# author: Joseph Lucero
# created on: 3 May 2018 08:23:46
# purpose: loading 3ddose files for plotting

from __future__ import division

from numpy import where, logical_and

from gzip import open as gOpen

import numpy


def position_to_index(p, A):
    """
    Description:
    Checks entire position array A until it finds the 
    Inputs:

    Outputs:
    """
    if p <= 0:
        return where(
            logical_and(p >= A[:-1], p < A[1:])
        )[0][0] + 1
    else:
        return where(
            logical_and(p >= A[:-1], p < A[1:])
        )[0][0]

class DoseFile(object):

    # copied from github/christopherpoole/3DDose

    def __init__(self, file_name, load_uncertainty=False):
        """
        Attempts to detect the dose file etension automatically. If an unknown
        extension is detected, loads a .3ddose file by default.
        """
        if file_name[-3:] == 'npz':
            self._load_npz(file_name)
        else:
            self._load_3ddose(file_name, load_uncertainty)
    
    def _load_npz(self, file_name):
        data = numpy.load(file_name)
        self.dose = data['dose']
        self.uncertainty = data['uncertainty']
        self.positions = [data['x_positions'], data['y_positions'], data['z_positions']]
        self.spacing = [numpy.diff(p) for p in self.positions]
        self.resolution = [s[0] for s in self.spacing if s.all()] 
        
        self.shape = self.dose.shape
        self.size = self.dose.size
    
    def _load_3ddose(self, file_name, load_uncertainty=False):
        if file_name[-3:] == '.gz':
            with gOpen(file_name, 'rb') as data_file:
                data = data_file.read().split('\n')
        else:
            data = file(file_name).read().split('\n')
        
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
        positions = [positions[2], positions[1], positions[0]]
        
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
        self.dose = self.dose.reshape((self.shape)).transpose()
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
                ).transpose()
            assert self.uncertainty.size == self.size, \
            "uncertainty array size does not match that specified."
            
    def dump(self, file_name):
        numpy.savez(file_name, dose=self.dose, uncertainty=self.uncertainty,
            x_positions=self.positions[0], y_positions=self.positions[1],
            z_positions=self.positions[2])

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

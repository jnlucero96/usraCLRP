#!/usr/bin/env
# author: Joseph Lucero
# created on: 12 Jun 2018 08:23:46
# purpose: calculate conversion factors needed to convert egs_brachy output to 
# desired output


from __future__ import division


def get_conversion_factor(air_kerma_strength, air_kerma_per_hist, max_indiv_dwell):
    """
    Description:
    Takes max air kerma strength, air kerma strength per history and max 
    individual dwell time to calculate the conversion factor needed to turn 
    dose per history output of egs_brachy to absolute dose in (Gy)

    Inputs:
    :param air_kerma_strength: Experimentally measured air kerma strength 
                               in units of (Gy cm^{2} h^{-1}) 
    :type air_kerma_strength: float
    :param air_kerma_per_hist: Simulation measured air kerma per history 
                               in units of (Gy cm^{2} history^{-1})
    :type air_kerma_per_hist: float
    :param max_indiv_dwell: Largest dwell time in treatment in units of (hours)
    :type max_indiv_dwell: float

    Outputs:
    :param conversion_factor: Conversion factor that turns dose per history to 
                              absolute dose 
    :type conversion_factor: float
    """

    return (air_kerma_strength * max_indiv_dwell) / air_kerma_per_hist

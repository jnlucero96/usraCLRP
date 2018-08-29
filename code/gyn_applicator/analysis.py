#!/usr/bin/env python2
# author: Joseph Lucero
# created on: 3 May 2018 08:23:46
# purpose: plotting 3ddose files

from __future__ import division

from sys import argv, exit
from os import getcwd

from numpy import (
    linspace, zeros_like, histogram, arange, sqrt, isnan, array, meshgrid,
    empty, nan_to_num, loadtxt
    )
from scipy.interpolate import RegularGridInterpolator as RGI
from py3ddose import DoseFile, position_to_index
from normalize import get_conversion_factor

from matplotlib.cm import get_cmap, tab10
from matplotlib.lines import Line2D
from matplotlib.style import use
use('seaborn-paper')
from matplotlib.pyplot import subplots, close, figure
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

def generate_tdvh_mlwa():

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' # for home
    target_dir = '/Users/student/research/results_not_to_git' # for work

    file_no_shield_0point5mm = target_dir + \
        '/mlwa_0shield_0.5mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_no_shield_1point0mm = target_dir + \
        '/mlwa_0shield_1mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_no_shield_2point0mm = target_dir + \
        '/mlwa_0shield_2mm_sim.phantom_wo_applicator_wo_box.3ddose'

    file_90_shield_0point5mm = target_dir + \
        '/mlwa_90shield_0.5mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_90_shield_1point0mm = target_dir + \
        '/mlwa_90shield_1mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_90_shield_2point0mm = target_dir + \
        '/mlwa_90shield_2mm_sim.phantom_wo_applicator_wo_box.3ddose'

    file_180_shield_0point5mm = target_dir + \
        '/mlwa_180shield_0.5mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_180_shield_1point0mm = target_dir + \
        '/mlwa_180shield_1mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_180_shield_2point0mm = target_dir + \
        '/mlwa_180shield_2mm_sim.phantom_wo_applicator_wo_box.3ddose'

    file_270_shield_0point5mm = target_dir + \
        '/mlwa_270shield_0.5mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_270_shield_1point0mm = target_dir + \
        '/mlwa_270shield_1mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_270_shield_2point0mm = target_dir + \
        '/mlwa_270shield_2mm_sim.phantom_wo_applicator_wo_box.3ddose'

    file_dict = {
        (0, 0): file_no_shield_0point5mm,
        (0, 1): file_90_shield_0point5mm,
        (0, 2): file_180_shield_0point5mm,
        (0, 3): file_270_shield_0point5mm,

        (1, 0): file_no_shield_1point0mm,
        (1, 1): file_90_shield_1point0mm,
        (1, 2): file_180_shield_1point0mm,
        (1, 3): file_270_shield_1point0mm,

        (2, 0): file_no_shield_2point0mm,
        (2, 1): file_90_shield_2point0mm,
        (2, 2): file_180_shield_2point0mm,
        (2, 3): file_270_shield_2point0mm
    }

    color_list = [
        'darkgreen',
        'purple',
        'darkorange',
        'saddlebrown'
    ]

    label_list = [
        'Unshielded',
        '90 degree shield',
        '180 degree shield',
        '270 degree shield'
    ]

    title_list = [
        'Voxel size = (0.5 mm)^{3}',
        'Voxel size = (1.0 mm)^{3}',
        'Voxel size = (2.0 mm)^{3}'
    ]

    fig, ax = subplots(3, 1, figsize=(10, 10), sharex='all', sharey='all')
    fig2, ax2 = subplots(3, 1, figsize=(10, 10), sharex='all', sharey='all')

    for index1 in xrange(3):
        for index2 in xrange(4):

            main_data = DoseFile(file_dict[(index1,index2)])

            # scale to absolute dose using maximum individual dwell time
            data_flat = main_data.dose.flatten() * 8.2573429808917e13 
            
            n, bins, __ = ax[index1].hist(
                data_flat,
                bins=500,
                color=color_list[index1],
                label=label_list[index1],
                weights=zeros_like(data_flat) + 1. / data_flat.size * 100,
                alpha=0.4
                )

            n_cum_base = n[::-1].cumsum()[::-1]

            ax2[index1].loglog(
                bins[:-1], n_cum_base,
                color=color_list[index1],
                label=label_list[index1]
            )
        
    for i in xrange(3):
        ax[i].legend(loc=0,prop={'size':10})
        ax[i].set_title(title_list[i],fontsize=20)
        ax[i].xaxis.set_tick_params(labelsize=17)
        ax[i].yaxis.set_tick_params(labelsize=17)
        
        ax2[i].legend(loc=0,prop={'size':10})
        ax2[i].set_title(title_list[i],fontsize=20)
        ax2[i].xaxis.set_tick_params(labelsize=17)
        ax2[i].yaxis.set_tick_params(labelsize=17)

    fig.text(
        0.03, 0.51, 'Volume (%)',
        fontsize=22, rotation='vertical', va='center'
    )
    fig.text(
        0.46, 0.03, 'Dose (Gy)', fontsize=22, va='center'
    )
    fig.text(
        0.52, 0.95,
        'Dose Volume Histograms \n With Volume Correction; ncase = 1E9',
        fontsize=22, va='center', ha='center'
    )

    fig.tight_layout()

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.88     # the top of the subplots of the figure
    # wspace = 0.2  # the amount of width for blank space between subplots
    # hspace = 0.2  # the amount of height for white space between subplotss

    fig.subplots_adjust(left=left, top=top, right=right, bottom=bottom)
    fig.savefig(pwd + '/dose_volume_histogram.pdf')
    
    fig2.text(
        0.03, 0.51, 'Volume (%)',
        fontsize=22, rotation='vertical', va='center'
    )
    fig2.text(
        0.45, 0.03, 'Dose (Gy)', fontsize=22, va='center'
    )
    fig2.text(
        0.52, 0.95,
        'Dose Volume Histograms \n With Volume Correction; ncase = 1E9',
        fontsize=22, va='center', ha='center'
    )

    fig2.tight_layout()

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.88     # the top of the subplots of the figure
    # wspace = 0.2  # the amount of width for blank space between subplots
    # hspace = 0.2  # the amount of height for white space between subplotss

    fig2.subplots_adjust(left=left, top=top, right=right, bottom=bottom)
    fig2.savefig(pwd + '/cumulative_dose_volume_histogram.pdf')

def dose_position_plots():
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots 
    from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' for home
    target_dir = '/Users/student/research/results_not_to_git' # for work

    file_list = [
        target_dir +
        '/mlwa_30mmOut_0shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        target_dir +
        '/mlwa_30mmOut_90shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        target_dir +
        '/mlwa_30mmOut_180shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        target_dir +
        '/mlwa_30mmOut_270shield_{0}_sim.phantom_wo_applicator.3ddose.gz'
    ]

    air_kerma_strength = 326.05715
    air_kerma_per_history = 1.1584E-13
    max_dwell_time = 0.02917

    vox_size_list = [
        # '0.5mm',
        '1pt0mm',
        '2pt0mm'
    ]

    shield_type_lst = [
        '0 shield',
        '90 shield',
        '180 shield',
        '270 shield'
    ]

    title_list = [
        'X - axis',
        'Y - axis',
        'Z - axis'
    ]

    for fig_index in xrange(len(vox_size_list)):

        fig, ax = subplots(
            3, len(file_list), figsize=(10, 10),
            sharex='all', sharey='all'
        )

        for index2 in xrange(len(file_list)):  # iterate through shield types

            full_data = DoseFile(
                file_list[index2].format(vox_size_list[fig_index]),
                load_uncertainty=True
            )

            Nx, Ny, Nz = full_data.shape

            # scale to absolute dose 
            full_data.dose *= get_conversion_factor(air_kerma_strength,
                air_kerma_per_history, max_dwell_time
            )
            full_data.uncertainty *= full_data.dose.__abs__()

            x_min, x_max = full_data.x_extent
            y_min, y_max = full_data.y_extent
            z_min, z_max = full_data.z_extent

            x_pos = array(full_data.positions[0])
            y_pos = array(full_data.positions[1])
            z_pos = array(full_data.positions[2])

            x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2.0
            y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2.0
            z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2.0

            for index1 in xrange(3):  # iterate through axes x, y, and z

                if index1 == 0:
                    ax[index1, index2].errorbar(
                        x_pos_mid, 
                        full_data.dose[
                            :, 
                            position_to_index(0.0, y_pos_mid),
                            position_to_index(0.0, z_pos_mid)
                            ],
                        yerr=full_data.uncertainty[
                            :, 
                            position_to_index(0.0, y_pos_mid), 
                            position_to_index(0.0, z_pos_mid)
                            ],
                        lw=3.0 
                    )
                elif index1 == 1:
                    ax[index1, index2].errorbar(
                        y_pos_mid,
                        full_data.dose[
                            position_to_index(0.0, x_pos_mid),
                            :,
                            position_to_index(0.0, z_pos_mid)
                            ],
                        yerr=full_data.uncertainty[
                            position_to_index(0.0, x_pos_mid), 
                            :, 
                            position_to_index(0.0, z_pos_mid)
                            ],
                        lw=3.0
                    )
                else:
                    ax[index1, index2].errorbar(
                        z_pos_mid,
                        full_data.dose[
                            position_to_index(0.0, x_pos_mid), 
                            position_to_index(0.0, y_pos_mid),
                            :
                            ],
                        yerr=full_data.uncertainty[
                            position_to_index(0.0, x_pos_mid), 
                            position_to_index(0.0, y_pos_mid), 
                            :
                            ],
                        lw=3.0
                    )

        for n in xrange(3):
            for m in xrange(len(file_list)):
                # ax[n, m].legend(loc=0,prop={'size':8})
                ax[n, m].grid(True)
                if n == 0:
                    ax[n, m].set_title(shield_type_lst[m])
                if m == len(file_list) - 1:
                    ax[n, m].set_ylabel(title_list[n],fontsize=20)
                    ax[n, m].yaxis.set_label_position('right')
                ax[n, m].xaxis.set_tick_params(labelsize=14)
                ax[n, m].yaxis.set_tick_params(labelsize=14)
                if n == 2:
                    ax[n, m].set_xticks(arange(-10, 10 + 1,5))
                    ax[n, m].set_xlim([-10,10])

        fig.text(
            0.01, 0.51, 'Dose (Gy)',
            fontsize=27, rotation='vertical', va='center', ha='center'
        )
        fig.text(
            0.51, 0.03, 'Position (cm)', fontsize=27, va='center', ha='center'
        )
        fig.text(
            0.52, 0.95,
            'Absolute Dose vs. Position \n With Volume Correction; ncase = 1E9',
            fontsize=27, va='center', ha='center'
        )
        fig.tight_layout()

        left = 0.125  # the left side of the subplots of the figure
        right = 0.95    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.87     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width for blank space between subplots
        # hspace = 0.2  # the amount of height for white space between subplots

        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig.savefig(
            pwd + '/dosage_comparison_' + vox_size_list[fig_index] + '.pdf'
            )   

def rel_dose_position_plot():
    """
    Description:
    Takes any number of .3ddose files and plots a figure similar to Fig. 5 of
    Lymperopoulou et al. (2004)

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' for home
    target_dir = '/Users/student/research/results_not_to_git'  # for work

    outer_diams = [
        # '25', 
        '30', 
        # '35', 
        # '40'
    ]

    diameter_to_radius = {
        # '25': 1.25,
        '30': 1.50,
        # '35': 1.75,
        # '40': 2.00
    }

    voxel_type = [
        '1pt0',
        '2pt0'
    ]

    for diameter in outer_diams:
        for vox_type in voxel_type:

            unshielded_file = target_dir + \
                '/mlwa_{0}mmOut_0shield_{1}mm'.format(diameter, vox_type) + \
                '_sim.phantom_wo_applicator.3ddose.gz'

            shielded_file_90 = target_dir + \
                '/mlwa_{0}mmOut_90shield_{1}mm'.format(diameter, vox_type) + \
                '_sim.phantom_wo_applicator.3ddose.gz'

            shielded_file_180 = target_dir + \
                '/mlwa_{0}mmOut_180shield_{1}mm'.format(diameter, vox_type) + \
                '_sim.phantom_wo_applicator.3ddose.gz'

            shielded_file_270 = target_dir + \
                '/mlwa_{0}mmOut_270shield_{1}mm'.format(diameter, vox_type) + \
                '_sim.phantom_wo_applicator.3ddose.gz'

            ref_lymper_90 = loadtxt(
                '/Users/student/research/reference/lymperopoulou_reference_data'
                + '/rel_dose_plot_compiled_90.dat'
            )
            ref_lymper_180 = loadtxt(
                '/Users/student/research/reference/lymperopoulou_reference_data'
                + '/rel_dose_plot_compiled_180.dat'
            )
            ref_lymper_270 = loadtxt(
                '/Users/student/research/reference/lymperopoulou_reference_data'
                + '/rel_dose_plot_compiled_270.dat'
            )

            if vox_type == '1pt0':
                stride = 4
            else:
                stride = 2            

            fig, ax = subplots(
                1, 1, figsize=(10, 7)
                )
            ax_twin = ax.twinx()

            ax.scatter(
                ref_lymper_90[:, 0], ref_lymper_90[:, 1], color='darkgreen',
                marker='o', s=14, 
                label=r'$90^{\circ}$ Ref.'
                )
            ax_twin.scatter(
                ref_lymper_90[:, 0], ref_lymper_90[:, 1], color='darkgreen',
                marker='o', s=14
                )

            ax.scatter(
                ref_lymper_180[:, 0], ref_lymper_180[:, 1], color='purple',
                marker='s', s=14, 
                label=r'$180^{\circ}$ Ref.'
                )
            ax_twin.scatter(
                ref_lymper_180[:, 0], ref_lymper_180[:, 1], color='purple',
                marker='s', s=14
                )

            ax.scatter(
                ref_lymper_270[:, 0], ref_lymper_270[:, 1], color='darkorange',
                marker='^', s=14, 
                label=r'$180^{\circ}$ Ref.'
                )
            ax_twin.scatter(
                ref_lymper_270[:, 0], ref_lymper_270[:, 1], color='darkorange',
                marker='^', s=14
                )

            unshielded_full_data = DoseFile(unshielded_file, load_uncertainty=True)
            shield_90_data = DoseFile(shielded_file_90, load_uncertainty=True)
            shield_180_data = DoseFile(shielded_file_180, load_uncertainty=True)
            shield_270_data = DoseFile(shielded_file_270, load_uncertainty=True)

            # print unshielded_full_data.shape
            # print shield_90_data.shape
            # print shield_180_data.shape
            # print shield_270_data.shape
            # exit(0)

            x_min, x_max = unshielded_full_data.x_extent
            z_min, z_max = unshielded_full_data.z_extent

            x_pos = array(unshielded_full_data.positions[0])
            y_pos = array(unshielded_full_data.positions[1])
            z_pos = array(unshielded_full_data.positions[2])

            x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2.0
            y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2.0
            z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2.0

            error_90 = sqrt(
                (shield_90_data.uncertainty) ** 2
                + (unshielded_full_data.uncertainty) ** 2
            ) 
            error_180 = sqrt(
                (shield_180_data.uncertainty) ** 2
                + (unshielded_full_data.uncertainty) ** 2
            ) 
            error_270 = sqrt(
                (shield_270_data.uncertainty) ** 2
                + (unshielded_full_data.uncertainty) ** 2
            ) 

            unshielded_interpolated_dose_matrix = unshielded_full_data.dose
            shield_90_interpolated_dose_matrix = shield_90_data.dose
            shield_180_interpolated_dose_matrix = shield_180_data.dose
            shield_270_interpolated_dose_matrix = shield_270_data.dose

            shield_90_interpolated_dose_matrix /= unshielded_interpolated_dose_matrix
            shield_180_interpolated_dose_matrix /= unshielded_interpolated_dose_matrix
            shield_270_interpolated_dose_matrix /= unshielded_interpolated_dose_matrix

            error_90 *= shield_90_data.dose.__abs__()
            error_180 *= shield_180_data.dose.__abs__()
            error_270 *= shield_270_data.dose.__abs__()

            ax.errorbar(
                x_pos_mid[::stride],
                # shield_90_data.dose[
                shield_90_interpolated_dose_matrix[
                # error_90[
                    :, 
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                yerr=error_90[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                lw=2.0, label=r'$90^{\circ}$ MC', color='darkgreen',
                # markeredgecolor='darkgreen', marker='o',
                markeredgewidth=1,# markersize=10, markerfacecolor='None',
                elinewidth=1.5, capsize=2.5
            )
            ax_twin.errorbar(
                x_pos_mid[::stride],
                shield_90_interpolated_dose_matrix[
                    :, 
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                yerr=error_90[
                    :, 
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                lw=2.0, label=r'$90^{\circ}$', color='darkgreen',
                # markeredgecolor='darkgreen', marker='o',
                markeredgewidth=1,# markersize=10, markerfacecolor='None',
                elinewidth=1.5, capsize=2.5
            )

            ax.errorbar(
                x_pos_mid[::stride],
                shield_180_interpolated_dose_matrix[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                yerr=error_180[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                lw=2.0, label=r'$180^{\circ}$ MC', color='purple',
                # markeredgecolor='purple', marker='s',
                markeredgewidth=1,# markersize=10, markerfacecolor='None',
                elinewidth=1.5, capsize=2.5
            )
            ax_twin.errorbar(
                x_pos_mid[::stride],
                shield_180_interpolated_dose_matrix[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                yerr=error_180[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                lw=2.0, label=r'$180^{\circ}$', color='purple',
                # markeredgecolor='purple', marker='s',
                markeredgewidth=1,# markersize=10, markerfacecolor='None',
                elinewidth=1.5, capsize=2.5
            )

            ax.errorbar(
                x_pos_mid[::stride],
                shield_270_interpolated_dose_matrix[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                yerr=error_270[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                lw=2.0, label=r'$270^{\circ}$ MC', color='darkorange',
                # markeredgecolor='darkorange', marker='^',
                markeredgewidth=1,# markersize=10, markerfacecolor='None',
                elinewidth=1.5, capsize=2.5
            )
            ax_twin.errorbar(
                x_pos_mid[::stride],
                shield_270_interpolated_dose_matrix[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                yerr=error_270[
                    :,
                    position_to_index(0.0, y_pos_mid),
                    position_to_index(0.0, z_pos_mid)
                ][::stride],
                lw=2.0, label=r'$270^{\circ}$', color='darkorange',
                # markeredgecolor='darkorange', marker='^',
                markeredgewidth=1,# markersize=10, markerfacecolor='None',
                elinewidth=1.5, capsize=2.5
            )

            ax.legend(loc='upper left', prop={'size': 12}, ncol=2)
            ax.grid(False)
            ax.set_yticks(arange(0.80, 1.11, 0.05))
            ax.set_ylim([0.80, 1.10])
            ax.xaxis.set_tick_params(labelsize=14)
            ax.yaxis.set_tick_params(labelsize=14)
            ax.set_xticks(arange(-10, 10 + 1, 2))
            ax.set_xlim([-10, 10])
            ax.axhline(1.0, 0, 0.5, lw=1.5, color='gray')
            
            radius = diameter_to_radius[diameter]
            ax_twin.grid(False)
            ax_twin.set_ylim([0.0, 0.3])
            ax_twin.set_yticks(arange(0.0, 0.31, 0.05))
            ax_twin.yaxis.set_tick_params(labelsize=14)
            ax_twin.vlines([-radius,radius],0,1.0)
            ax_twin.fill_between([-radius, radius], 1.05, facecolor='lightgray')

            fig.text(
                0.02, 0.52, 'Relative Dose',
                fontsize=27, rotation='vertical', va='center', ha='center'
            )
            fig.text(
                0.98, 0.52, 'Effective Transmission',
                fontsize=27, rotation='vertical', va='center', ha='center'
            )
            fig.text(
                0.5, 0.51, 'A p p l i c a t o r',
                fontsize=35, rotation='vertical', va='center', ha='center'
            )
            fig.text(
                0.28, 0.025, 'x-axis (cm)', fontsize=20, va='center', ha='center'
            )
            fig.text(
                0.72, 0.025, 'x-axis (cm)', fontsize=20, va='center', ha='center'
            )
            fig.tight_layout()

            left = 0.11  # the left side of the subplots of the figure
            right = 0.89    # the right side of the subplots of the figure
            bottom = 0.09   # the bottom of the subplots of the figure
            top = 0.97     # the top of the subplots of the figure
            # wspace = 0.2  # the amount of width for blank space between subplots
            # hspace = 0.2  # the amount of height for white space between subplotss

            fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

            fig.savefig(
                pwd 
                + '/relative_dosage_comparison_'
                + '{0}mmOut_{1}mmVoxel.pdf'.format(diameter,vox_type)
                )

            close(fig)

def dose_inv_position_plots(interpolate=False, plot=False):
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots 
    from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' for home
    target_dir = '/Users/student/research/results_not_to_git'  # for work

    file_list = [
        # target_dir +
        # '/mlwa_0shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        # target_dir +
        # '/mlwa_90shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        # target_dir +
        # '/mlwa_180shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        # target_dir +
        # '/mlwa_270shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_25mmOut_{0}shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        # target_dir +
        # '/mlwa_30mmOut_{0}shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        # target_dir +
        # '/mlwa_30mmOut_shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        # target_dir +
        # '/mlwa_30mmOut_270shield_{0}_sim.phantom_wo_applicator.3ddose.gz'
        # target_dir +
        # '/tg43_{0}_sim.phantom_wo_box.3ddose',
        # target_dir +
        # '/tg43_applicator_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
    ]

    shield_type_lst = [
        '0'
        # '90',
        # '180',
        # '270'
        # '180'
    ]

    vox_size_lst = [
        # '0.5mm',
        '1mm',
        '2mm'
    ]
    vox_size_txt_lst = [
        # '0.5mm',
        '1pt0',
        '2pt0'
    ]

    title_list = [
        'y-cuts (|x| = {0}, z = {1})',
        'z-cuts (|x| = {0}, y = {1})'
    ]

    # Lymperopoulou
    # air_kerma_true = 326.05715
    # air_kerma_per_hist = 1.1584e-13
    # max_dwell_time = 0.02917

    # Ali-Cygler
    air_kerma_str = 431.6
    air_kerma_per_hist = 1.1517e-13
    max_dwell_time = 0.01352278

    for x_pos_desired in [1.51, 1.61]:
        print '=' * 40
        print "x-position =", x_pos_desired
        for y_pos_desired in [0.01]:
            print "y-position =", y_pos_desired
            for fig_index, vox_size in enumerate(vox_size_txt_lst):
                # iterate through shield types
                for index2, shield_type in enumerate(shield_type_lst):

                    full_data = DoseFile(
                        target_dir
                        + '/mlwa_25mmOut_'
                        + '{0}shield_{1}mm'.format(shield_type, vox_size)
                        + '_sim.phantom_wo_applicator.3ddose.gz',
                        load_uncertainty=True
                    )

                    # scale to maximum individual dwell time
                    full_data.dose *= get_conversion_factor(
                        air_kerma_str, air_kerma_per_hist, max_dwell_time
                    )

                    Nx, Ny, Nz = full_data.shape

                    x_min, x_max = full_data.x_extent
                    y_min, y_max = full_data.y_extent
                    z_min, z_max = full_data.z_extent

                    x_pos = array(full_data.positions[0])
                    y_pos = array(full_data.positions[1])
                    z_pos = array(full_data.positions[2])

                    x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2.0
                    y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2.0
                    z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2.0

                    if interpolate:

                        interpolated_dose = RGI(
                            (x_pos_mid, y_pos_mid, z_pos_mid),
                            full_data.dose,
                            bounds_error=False,
                            fill_value=None
                        )

                        points = meshgrid(x_pos, y_pos, z_pos)
                        flat = array([m.flatten() for m in points])

                        # something happens during the reshape that messes up the
                        # (x, y, z) ordering not quite sure why.
                        # Need the transpose in order to get proper (x, y, z)
                        # ordering
                        dose_matrix = interpolated_dose(
                            flat.transpose()
                        ).reshape(*points[0].shape).transpose((1, 0, 2))

                        z_depths = (
                                dose_matrix[
                                    position_to_index(x_pos_desired, x_pos),
                                    position_to_index(y_pos_desired, y_pos),
                                    :
                                ] / dose_matrix[
                                    position_to_index(-x_pos_desired, x_pos),
                                    position_to_index(y_pos_desired, y_pos),
                                    :
                                ]
                            )

                    else:

                        dose_matrix = full_data.dose
                        err_matrix = full_data.uncertainty * dose_matrix.__abs__()

                        z_depths = (
                            dose_matrix[
                                position_to_index(x_pos_desired, x_pos),
                                position_to_index(y_pos_desired, y_pos),
                                :
                            ] / dose_matrix[
                                position_to_index(-x_pos_desired,x_pos),
                                position_to_index(y_pos_desired, y_pos),
                                :
                            ]
                        )
                        z_depths_err = sqrt(
                            (err_matrix[
                                position_to_index(x_pos_desired, x_pos),
                                position_to_index(y_pos_desired, y_pos),
                                :
                            ] ** 2) * (err_matrix[
                                position_to_index(-x_pos_desired,x_pos),
                                position_to_index(y_pos_desired, y_pos),
                                :
                            ] ** 2)
                        ) * z_depths.__abs__()

                    print "For shield type:", shield_type, ", voxel size:", vox_size
                    print "mean =", z_depths[
                        position_to_index(-2.5, z_pos):position_to_index(2.6, z_pos)
                        ].mean(), ", std. deviation =", z_depths[
                            position_to_index(-2.5, z_pos):position_to_index(2.6, z_pos)
                        ].std()
                    print "At origin:", z_depths[position_to_index(0.0, z_pos)], z_depths_err[position_to_index(0.0, z_pos)]
                    print; 

                    if plot:
                        
                        fig = figure(figsize=(15, 10))
                        minor_locator = MultipleLocator(0.5)
                        minor_locator2 = MultipleLocator(0.5)
                        gs = GridSpec(ncols=2, nrows=2)

                        ax1 = fig.add_subplot(gs[0, 0])
                        ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
                        ax3 = fig.add_subplot(gs[1, :])

                        ax = [ax1, ax2, ax3]

                        ax[0].errorbar(
                            z_pos_mid,
                            dose_matrix[
                                position_to_index(x_pos_desired, x_pos),
                                position_to_index(y_pos_desired, y_pos),
                                :
                            ],
                            yerr=err_matrix[
                                position_to_index(x_pos_desired, x_pos),
                                position_to_index(y_pos_desired, y_pos),
                                :
                                ],
                            lw=3.0, capsize=2.0, elinewidth=2.0
                        )
                        ax[1].errorbar(
                            z_pos_mid,
                            dose_matrix[
                                position_to_index(-x_pos_desired, x_pos),
                                position_to_index(y_pos_desired, y_pos),
                                :
                            ],
                            yerr=err_matrix[
                                position_to_index(-x_pos_desired, x_pos),
                                position_to_index(y_pos_desired, y_pos),
                                :
                                ],
                            lw=3.0, capsize=2.0, elinewidth=2.0
                        )
                        ax[2].errorbar(
                            z_pos_mid,
                            z_depths * 100,
                            yerr=z_depths_err * 100,
                            lw=3.0, capsize=2.0, elinewidth=2.0
                        )

                        ax[0].set_ylabel(
                            r"$D_{\mathrm{shielded}}$ (Gy)", fontsize=17
                        )
                        ax[1].set_ylabel(
                            r"$D_{\mathrm{unshielded}}$ (Gy)", fontsize=17
                        )
                        ax[2].set_ylabel(
                            r"$D_{\mathrm{shielded}}\ /\ D_{\mathrm{unshielded}}$ (%)",
                            fontsize=17
                        )

                        for m in xrange(3):
                            ax[m].grid(True, which='both')
                            ax[m].set_xticks(arange(-10, 10 + 1, 2))
                            ax[m].xaxis.set_minor_locator(minor_locator)
                            if m != 2:
                                ax[m].yaxis.set_minor_locator(minor_locator2)
                            ax[m].xaxis.set_tick_params(labelsize=12)
                            ax[m].yaxis.set_tick_params(labelsize=12)
                            ax[m].set_xlim([-10, 10])

                        # fig.text(
                        #     0.01, 0.51, 'Ratio at depths',
                        #     fontsize=27, rotation='vertical', va='center'
                        # )
                        fig.text(
                            0.55, 0.03, 'Depth of cut (cm)', fontsize=27, va='center',
                            ha='center'
                        )
                        # fig.text(
                        #     0.52, 0.95,
                        #     'Dose Ratios vs. Depth \n With Volume Correction; ncase = 1E9',
                        #     fontsize=27, va='center', ha='center'
                        # )
                        fig.tight_layout()

                        left = 0.09  # the left side of the subplots of the figure
                        right = 0.97    # the right side of the subplots of the figure
                        bottom = 0.09   # the bottom of the subplots of the figure
                        top = 0.95     # the top of the subplots of the figure
                        # wspace = 0.2  # the amount of width for blank space between subplots
                        # hspace = 0.2  # the amount of height for white space between subplots

                        fig.subplots_adjust(
                            left=left, bottom=bottom, right=right, top=top)

                        fig.savefig(
                            pwd + '/dosage_inv_comparison_' + vox_size
                            + '_shield' + shield_type
                            + '_x' + str(x_pos_desired) + '_y' + str(y_pos_desired)
                            + '.pdf'
                        )

                        close(fig)

def isodose_plot(mode='mlwa'):
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots 
    from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' # for home
    target_dir = '/Users/student/research/results_not_to_git'  # for work

    air_kerma_true = 326.05715
    air_kerma_per_hist = 1.1584e-13
    max_dwell_time = 0.02917
    # air_kerma_true = 431.6
    # air_kerma_per_hist = 1.1517e-13
    # max_dwell_time = 0.01352278

    if mode == 'mlwa':
        target_file = '/mlwa_{1}mmOut_{0}shield_{2}pt0mm_sim.phantom_wo_' + \
            'applicator.3ddose.gz'
        base_title_str = 'Nucletron {0}mm diameter applicator \n Seed' + \
            '-as-source; ncase = 5E9'
    elif mode == 'tg43appl':
        target_file = '/tg43appl_{1}mmOut_{2}pt0mm_sim.phantom_wo_applicator_wo_box.3ddose.gz'
        base_title_str = 'Ghost {0}mm diameter applicator \n Seed' + \
            '-as-source; ncase = 5E9'
    elif mode == 'tg43pure':
        target_file = '/tg43_{2}pt0mm_sim.phantom_wo_box.3ddose.gz'
        base_title_str = 'Pure TG43 \n Seed-as-source; ncase = 5E9'
    else: 
        print "Mode not understood. Exiting now"
        exit(1)

    
    outer_diams = [
        '25'
        # '30'
        # '35',
        # '40'
    ]

    shield_type_lst = [
        '0',
        '90',
        '180',
        '270',
    ]

    vox_size_lst = [
        '1',
        '2',
    ]

    for diam_index, diameter in enumerate(outer_diams):

        for vox_index, vox_size in enumerate(vox_size_lst):

            title_str = base_title_str.format(diameter)

            fig, ax = subplots(
                2, 2, figsize=(12, 12),
                sharex='all', sharey='all'
            )
            fig2, ax2 = subplots(
                2, 2, figsize=(12, 12),
                sharex='all', sharey='all'
            )

            for shield_index, shield_type in enumerate(shield_type_lst):

                full_data = DoseFile(
                    target_dir + target_file.format(
                        shield_type, diameter, vox_size
                        )
                    )

                x_pos = array(full_data.positions[0])
                y_pos = array(full_data.positions[1])
                z_pos = array(full_data.positions[2])

                x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2.0
                y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2.0
                z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2.0

                if shield_index == 0:
                    ax_x, ax_y = 0, 0
                elif shield_index == 1:
                    ax_x, ax_y = 0, 1
                elif shield_index == 2:
                    ax_x, ax_y = 1, 0
                else:
                    ax_x, ax_y = 1, 1

                Nx, Ny, Nz = full_data.shape

                full_data.dose *= get_conversion_factor(
                    air_kerma_true, 
                    air_kerma_per_hist,
                    max_dwell_time

                )  # scale to maximum individual dwell time

                full_data.dose /= 5  # normalize to desired dose of 5 Gy
                full_data.dose *= 100  # express in percent. Should see 100% at x=-2cm

                xy_contour = ax[ax_x, ax_y].contourf(
                    x_pos_mid, y_pos_mid, 
                    # matplotlib plots column by row (instead of row by column)
                    # so transpose data array to account for this
                    full_data.dose[:, :, position_to_index(0.0, z_pos_mid)].transpose(),
                    arange(0, 110, 10),
                    # [5, 10, 20, 50, 100],
                    # cmap=get_cmap('gnuplot')
                )
                ax[ax_x, ax_y].set_title(
                    r'${0}$'.format(shield_type) + r'$^{\circ}$ shielding',
                    fontsize=20
                    )

                xz_contour = ax2[ax_x, ax_y].contourf(
                    x_pos_mid, z_pos_mid,
                    # matplotlib plots column by row (instead of row by column)
                    # so transpose data array to account for this
                    full_data.dose[:, position_to_index(0.0, y_pos_mid), :].transpose(),
                    arange(0, 110, 10),
                    # [5, 10, 20, 50, 100],
                    # cmap=get_cmap('gnuplot')
                )
                ax2[ax_x, ax_y].set_title(
                    r'${0}$'.format(shield_type) + r'$^{\circ}$ shielding',\
                    fontsize=20
                    )

            for n in xrange(2):
                for m in xrange(2):

                    # ax[n, m].grid(True)
                    ax[n, m].xaxis.set_tick_params(
                        labelsize=15, top=True, direction='in'
                        )
                    ax[n, m].yaxis.set_tick_params(
                        labelsize=15, right=True, direction='in'
                        )
                    ax[n, m].set_xticks(arange(-10, 10 + 1, 2))
                    ax[n, m].set_xlim([-10,10])
                    ax[n, m].set_yticks(arange(-10, 10 + 1, 2))
                    ax[n, m].set_ylim([-10, 10])
                    # ax[n, m].vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
                    # ax[n, m].hlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
                    # ax[n, m].vlines([0], -10, 10, linestyles='dashed',
                    #                 lw=2.0, color='darkgoldenrod')
                    # ax[n, m].hlines([0], -10, 10, linestyles='dashed',
                    #                 lw=2.0, color='darkgoldenrod')

                    # ax2[n, m].grid(True)
                    ax2[n, m].xaxis.set_tick_params(
                        labelsize=15, top=True, direction='in'
                        )
                    ax2[n, m].yaxis.set_tick_params(
                        labelsize=15, right=True, direction='in'
                        )
                    ax2[n, m].set_xticks(arange(-10, 10 + 1, 2))
                    ax2[n, m].set_xlim([-10, 10])
                    ax2[n, m].set_yticks(arange(-10, 10 + 1, 2))
                    ax2[n, m].set_ylim([-10, 10])

                    # ax2[n, m].vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
                    # ax2[n, m].vlines([0], -10, 10, linestyles='dashed', lw=2.0, color='darkgoldenrod')
                    # ax2[n, m].hlines([0], -10, 10, linestyles='dashed',
                    #                  lw=2.0, color='darkgoldenrod')
                    
            fig.text(
                0.01, 0.51, 'y-axis (cm)',
                fontsize=30, rotation='vertical', va='center'
            )
            fig.text(
                0.51, 0.03, 'x-axis (cm)', fontsize=30, va='center',ha='center'
            )

            cax = fig.add_axes([0.90, 0.09, 0.01, 0.86])
            cbar1 = fig.colorbar(
                xy_contour, cax=cax, orientation='vertical',
                ax=ax
                # ticks=[0, 20, 40, 60, 80, 100]
            )
            cbar1.set_label('Percentage Isodose (%)', fontsize=24)
            # cbar1.set_label('Dose (Gy)', fontsize=24)
            # cbar1.ax.set_yticklabels([0, 1, 2, 3, 4, 5])
            cbar1.ax.tick_params(labelsize=20)
            fig.tight_layout()

            left = 0.1  # the left side of the subplots of the figure
            right = 0.89    # the right side of the subplots of the figure
            bottom = 0.09   # the bottom of the subplots of the figure
            top = 0.95     # the top of the subplots of the figure
            # wspace = 0.2  # the amount of width reserved for blank space between subplots
            # hspace = 0.2  # the amount of height reserved for white space between subplotss

            fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

            fig.savefig(
                pwd + '/xy_isodose_profile_' + diameter + 'mmOut_' 
                + vox_size + 'mm.pdf'
            )

            fig2.text(
                0.01, 0.51, 'z-axis (cm)',
                fontsize=30, rotation='vertical', va='center'
            )
            fig2.text(
                0.43, 0.03, 'x-axis (cm)', fontsize=30, va='center'
            )
            cax2 = fig2.add_axes([0.90, 0.09, 0.01, 0.86])
            cbar2 = fig2.colorbar(
                xz_contour, cax=cax2, orientation='vertical',
                ax=ax
                # ticks=[0, 20, 40, 60, 80, 100]
            )
            cbar2.set_label('Percentage Isodose (%)', fontsize=20)
            # cbar2.set_label('Dose (Gy)', fontsize=20)
            # cbar2.ax.set_yticklabels([0, 1, 2, 3, 4, 5])
            cbar2.ax.tick_params(labelsize=20)

            fig2.tight_layout()

            left = 0.1  # the left side of the subplots of the figure
            right = 0.89    # the right side of the subplots of the figure
            bottom = 0.09   # the bottom of the subplots of the figure
            top = 0.95     # the top of the subplots of the figure
            # wspace = 0.2  # the amount of width reserved for blank space between subplots
            # hspace = 0.2  # the amount of height reserved for white space between subplotss

            fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

            fig2.savefig(
                pwd + '/xz_isodose_profile_' + diameter + 'mmOut_' 
                + vox_size + 'mm.pdf'
            )

def isodose_plot_compare():
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots 
    from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' # for home
    target_dir = '/Users/student/research/results_not_to_git'  # for work
    target_ref_dir = '/Users/student/research/reference/' \
        + 'lymperopoulou_reference_data'

    air_kerma_true = 326.05715
    air_kerma_per_hist = 1.1584e-13
    max_dwell_time = 0.02917
    # air_kerma_true = 431.6
    # air_kerma_per_hist = 1.1517e-13
    # max_dwell_time = 0.01352278

    target_file = '/mlwa_{1}mmOut_{0}shield_{2}pt0mm_sim.phantom_wo_' + \
        'applicator.3ddose.gz'
    reference_file_xy = '/xy_isodose_contour_{0}shield_compiled.dat'
    reference_file_xz = '/xz_isodose_contour_{0}shield_compiled.dat'

    outer_diams = [
        # '25'
        '30'
        # '35',
        # '40'
    ]

    shield_type_lst = [
        '0',
        '90',
        '180',
        '270',
    ]

    vox_size_lst = [
        '1'
        # '2'
    ]

    for diam_index, diameter in enumerate(outer_diams):

        for vox_index, vox_size in enumerate(vox_size_lst):

            fig, ax = subplots(
                2, 2, figsize=(12, 12),
                sharex='all', sharey='all'
            )
            fig2, ax2 = subplots(
                2, 2, figsize=(12, 12),
                sharex='all', sharey='all'
            )

            for shield_index, shield_type in enumerate(shield_type_lst):

                full_data = DoseFile(
                    target_dir + target_file.format(
                        shield_type, diameter, vox_size
                    )
                )

                x_pos = array(full_data.positions[0])
                y_pos = array(full_data.positions[1])
                z_pos = array(full_data.positions[2])

                x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2.0
                y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2.0
                z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2.0

                ax_x = shield_index // 2
                ax_y = shield_index % 2

                Nx, Ny, Nz = full_data.shape

                full_data.dose *= get_conversion_factor(
                    air_kerma_true,
                    air_kerma_per_hist,
                    max_dwell_time

                )  # scale to maximum individual dwell time

                full_data.dose /= 5  # normalize to desired dose of 5 Gy
                full_data.dose *= 100  # express in percent. Should see 100% at x=-2cm

                xy_contour = ax[ax_x, ax_y].contour(
                    x_pos_mid, y_pos_mid,
                    # matplotlib plots column by row (instead of row by column)
                    # so transpose data array to account for this
                    full_data.dose[:, :, position_to_index(0.0, z_pos)].transpose(),
                    [5, 10, 20, 50, 100],
                    colors=tab10(linspace(0, 1, 5))
                )
                ax[ax_x, ax_y].set_title(
                    r'${0}$'.format(shield_type) + r'$^{\circ}$ shielding',
                    fontsize=20
                )

                xz_contour = ax2[ax_x, ax_y].contour(
                    x_pos_mid, z_pos_mid,
                    # matplotlib plots column by row (instead of row by column)
                    # so transpose data array to account for this
                    full_data.dose[:, position_to_index(0.0, y_pos), :].transpose(),
                    [5, 10, 20, 50, 100],
                    colors=tab10(linspace(0, 1, 5))
                )
                ax2[ax_x, ax_y].set_title(
                    r'${0}$'.format(shield_type) + r'$^{\circ}$ shielding',
                    fontsize=20
                )

                if shield_type != '0':
                    ref_data_xy = loadtxt(
                        target_ref_dir + reference_file_xy.format(shield_type)
                    )
                    ref_data_xz = loadtxt(
                        target_ref_dir + reference_file_xz.format(shield_type)
                    )
                    ax[ax_x, ax_y].scatter(
                        ref_data_xy[:, 0], ref_data_xy[:, 1], marker='x', s=14,
                        color='black'
                    )
                    ax2[ax_x, ax_y].scatter(
                        ref_data_xz[:, 0], ref_data_xz[:, 1], marker='x', s=14,
                        color='black'
                    )

            for n in xrange(2):
                for m in xrange(2):

                    ax[n, m].xaxis.set_tick_params(
                        labelsize=15, top=True, direction='in'
                    )
                    ax[n, m].yaxis.set_tick_params(
                        labelsize=15, right=True, direction='in'
                    )
                    ax[n, m].set_xticks(arange(-10, 10 + 1, 2))
                    ax[n, m].set_xlim([-10, 10])
                    ax[n, m].set_yticks(arange(-10, 10 + 1, 2))
                    ax[n, m].set_ylim([-10, 10])

                    ax2[n, m].xaxis.set_tick_params(
                        labelsize=15, top=True, direction='in'
                    )
                    ax2[n, m].yaxis.set_tick_params(
                        labelsize=15, right=True, direction='in'
                    )
                    ax2[n, m].set_xticks(arange(-10, 10 + 1, 2))
                    ax2[n, m].set_xlim([-10, 10])
                    ax2[n, m].set_yticks(arange(-10, 10 + 1, 2))
                    ax2[n, m].set_ylim([-10, 10])

            fig.text(
                0.025, 0.51, 'y-axis (cm)', fontsize=30, rotation='vertical', 
                va='center', ha='center'
            )
            fig.text(
                0.54, 0.03, 'x-axis (cm)', fontsize=30, 
                va='center', ha='center'
            )
            fig.text(
                0.54, 0.97, 'Percent Isodose (%)', fontsize=30,
                va='center', ha='center'
            )

            fig.tight_layout()

            left = 0.1  # the left side of the subplots of the figure
            right = 0.95    # the right side of the subplots of the figure
            bottom = 0.09   # the bottom of the subplots of the figure
            top = 0.86     # the top of the subplots of the figure

            fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

            handles_list = [
                Line2D([], [], color=c, lw=3.0)
                for c in tab10(linspace(0, 1, 5))
            ]
            label_list = [
                '5', '10', '20', '50', '100'
            ]
            first_legend = fig.legend(
                handles_list, label_list, 
                loc='center', bbox_to_anchor=(0.53, 0.92),
                prop={'size': 22}, ncol=5,
                handlelength=2, frameon=False
            )
            fig.gca().add_artist(first_legend)

            fig.savefig(
                pwd + '/xy_isodose_profile_' + diameter + 'mmOut_'
                + vox_size + 'mm_w_ref.pdf'
            )

            fig2.text(
                0.025, 0.51, 'z-axis (cm)',fontsize=30, rotation='vertical', 
                va='center', ha='center'
            )
            fig2.text(
                0.54, 0.03, 'x-axis (cm)', fontsize=30, 
                va='center', ha='center'
            )
            fig2.text(
                0.54, 0.97, 'Percent Isodose (%)', fontsize=30,
                va='center', ha='center'
            )

            fig2.tight_layout()

            left = 0.1  # the left side of the subplots of the figure
            right = 0.95    # the right side of the subplots of the figure
            bottom = 0.09   # the bottom of the subplots of the figure
            top = 0.86     # the top of the subplots of the figure

            fig2.subplots_adjust(left=left, bottom=bottom,
                                 right=right, top=top)

            second_legend = fig2.legend(
                handles_list, label_list,
                loc='center', bbox_to_anchor=(0.53, 0.92),
                prop={'size': 22}, ncol=5,
                handlelength=2, frameon=False
            )
            fig2.gca().add_artist(second_legend)

            fig2.savefig(
                pwd + '/xz_isodose_profile_' + diameter + 'mmOut_'
                + vox_size + 'mm_w_ref.pdf'
            )

def tg43_mbdca_comparison_isodose_plot(explicit_contour=False):
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots 
    from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' # for home
    target_dir = '/Users/student/research/results_not_to_git'  # for work

    air_kerma_true = 326.05715
    air_kerma_per_hist = 1.1584e-13
    max_dwell_time = 0.02917

    # tg43_target_file = '/tg43_{0}pt0mm_sim.phantom_wo_box.3ddose.gz'  # for pure TG43 results
    # for ghost TG43 results
    tg43_target_file = '/tg43appl_{1}mmOut_{0}pt0mm_sim.phantom_wo_applicator_wo_box.3ddose.gz'
    mlwa_target_file = '/mlwa_{1}mmOut_{0}shield_{2}pt0mm_sim.phantom_wo_applicator.3ddose.gz'
    base_title_str = 'TG43-MBDCA comparison \n {0}mm applicator; ncase = 5E9'

    outer_diams = [
        # '25',
        '30'
        # '35',
        # '40'
    ]

    shield_type_lst = [
        '0',
        '90',
        '180',
        '270',
    ]

    vox_size_lst = [
        '1',
        '2',
    ]

    dose_scale_factor = get_conversion_factor(
        air_kerma_true,
        air_kerma_per_hist,
        max_dwell_time
    )  # scale to maximum individual dwell time

    for diam_index, diameter in enumerate(outer_diams):

        for vox_index, vox_size in enumerate(vox_size_lst):

            title_str = base_title_str.format(diameter)

            fig, ax = subplots(
                2, 2, figsize=(10, 10),
                sharex='all', sharey='all'
            )
            fig2, ax2 = subplots(
                2, 2, figsize=(10, 10),
                sharex='all', sharey='all'
            )

            tg43_full_data = DoseFile(
                target_dir + tg43_target_file.format(vox_size, diameter)
            )

            for shield_index, shield_type in enumerate(shield_type_lst):

                mlwa_full_data = DoseFile(
                    target_dir
                    + mlwa_target_file.format(shield_type, diameter, vox_size)
                )

                x_pos = array(tg43_full_data.positions[0])
                y_pos = array(tg43_full_data.positions[1])
                z_pos = array(tg43_full_data.positions[2])

                x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2.0
                y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2.0
                z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2.0

                ax_x = shield_index // 2
                ax_y = shield_index % 2

                tg43_dose = tg43_full_data.dose * dose_scale_factor
                mlwa_dose = mlwa_full_data.dose * dose_scale_factor

                # calculate percentage difference between MBDCA calculation and
                # tg43. Turn all NaNs to 0s and all inf to largest number that
                # is held by the double type
                per_diff = nan_to_num(
                    ((tg43_dose - mlwa_dose) / mlwa_dose) * 100
                    )

                xy_contour = ax[ax_x, ax_y].contourf(
                    x_pos_mid, y_pos_mid,
                    # matplotlib plots column by row (instead of row by column)
                    # so transpose data array to account for this
                    per_diff[:, :, position_to_index(
                        0.0, z_pos_mid)].transpose(),
                    arange(-150, 1200+1, 150),
                    # [5, 10, 20, 50, 100],
                    cmap=get_cmap('hot')
                )
                ax[ax_x, ax_y].set_title(
                    r'${0}$'.format(shield_type) + r'$^{\circ}$ shielding',
                    fontsize=15
                )

                xz_contour = ax2[ax_x, ax_y].contourf(
                    x_pos_mid, z_pos_mid,
                    # matplotlib plots column by row (instead of row by column)
                    # so transpose data array to account for this
                    per_diff[:, position_to_index(
                        0.0, y_pos_mid), :].transpose(),
                    arange(-150, 1200+1, 150),
                    # [5, 10, 20, 50, 100],
                    cmap=get_cmap('hot')
                )
                ax2[ax_x, ax_y].set_title(
                    r'${0}$'.format(shield_type) + r'$^{\circ}$ shielding',
                    fontsize=15
                )

                if explicit_contour:
                    explicit_XY = ax[ax_x, ax_y].contour(
                        x_pos_mid, y_pos_mid,
                        # matplotlib plots column by row (instead of row by column)
                        # so transpose data array to account for this
                        per_diff[:, :, position_to_index(
                            0.0, z_pos_mid)].transpose(),
                        # arange(0, 110, 10),
                        [10, 20, 50, 100],
                        cmap=get_cmap('Set1')
                    )
                    # ax[ax_x, ax_y].clabel(explicit_XY, inline=True, fontsize=4)
                    explicit_XZ = ax2[ax_x, ax_y].contour(
                        x_pos_mid, z_pos_mid,
                        # matplotlib plots column by row (instead of row by column)
                        # so transpose data array to account for this
                        per_diff[:, position_to_index(
                            0.0, y_pos_mid), :].transpose(),
                        # arange(0, 110, 10),
                        [10, 20, 50, 100],
                        cmap=get_cmap('Set1')
                    )
                    # ax2[ax_x, ax_y].clabel(explicit_XZ, inline=True, fontsize=4)

            for n in xrange(2):
                for m in xrange(2):

                    # ax[n, m].grid(True)
                    ax[n, m].xaxis.set_tick_params(
                        labelsize=15, top=True, direction='in'
                    )
                    ax[n, m].yaxis.set_tick_params(
                        labelsize=15, right=True, direction='in'
                    )
                    ax[n, m].set_xticks(arange(-10, 10 + 1, 2))
                    ax[n, m].set_xlim([-10, 10])
                    ax[n, m].set_yticks(arange(-10, 10 + 1, 2))
                    ax[n, m].set_ylim([-10, 10])
                    # ax[n, m].vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
                    # ax[n, m].hlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
                    # ax[n, m].vlines([0], -10, 10, linestyles='dashed',
                    #                 lw=2.0, color='darkgoldenrod')
                    # ax[n, m].hlines([0], -10, 10, linestyles='dashed',
                    #                 lw=2.0, color='darkgoldenrod')

                    # ax2[n, m].grid(True)
                    ax2[n, m].xaxis.set_tick_params(
                        labelsize=15, top=True, direction='in'
                    )
                    ax2[n, m].yaxis.set_tick_params(
                        labelsize=15, right=True, direction='in'
                    )
                    ax2[n, m].set_xticks(arange(-10, 10 + 1, 2))
                    ax2[n, m].set_xlim([-10, 10])
                    ax2[n, m].set_yticks(arange(-10, 10 + 1, 2))
                    ax2[n, m].set_ylim([-10, 10])

                    # ax2[n, m].vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
                    # ax2[n, m].vlines([0], -10, 10, linestyles='dashed', lw=2.0, color='darkgoldenrod')
                    # ax2[n, m].hlines([0], -10, 10, linestyles='dashed',
                    #                  lw=2.0, color='darkgoldenrod')

            fig.text(
                0.01, 0.51, 'y-axis (cm)',
                fontsize=27, rotation='vertical', va='center'
            )
            fig.text(
                0.43, 0.03, 'x-axis (cm)', fontsize=27, va='center'
            )

            cax = fig.add_axes([0.91, 0.09, 0.01, 0.86])
            cbar1 = fig.colorbar(
                xy_contour, cax=cax, orientation='vertical',
                ax=ax
            )
            cbar1.set_label('Percentage difference (%)', fontsize=17)
            cbar1.ax.tick_params(labelsize=10)
            fig.tight_layout()

            left = 0.1  # the left side of the subplots of the figure
            right = 0.90    # the right side of the subplots of the figure
            bottom = 0.09   # the bottom of the subplots of the figure
            top = 0.95     # the top of the subplots of the figure
            # wspace = 0.2  # the amount of width reserved for blank space between subplots
            # hspace = 0.2  # the amount of height reserved for white space between subplotss

            fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

            fig.savefig(
                pwd + '/tg43appl_mbdca_compare_xy_isodose_profile_' + diameter + 'mmOut_'
                + vox_size + 'mm.pdf'
            )

            fig2.text(
                0.01, 0.51, 'z-axis (cm)',
                fontsize=27, rotation='vertical', va='center'
            )
            fig2.text(
                0.43, 0.03, 'x-axis (cm)', fontsize=27, va='center'
            )
            cax2 = fig2.add_axes([0.91, 0.09, 0.01, 0.86])
            cbar2 = fig2.colorbar(
                xy_contour, cax=cax2, orientation='vertical',
                ax=ax
            )
            cbar2.set_label('Percentage difference (%)', fontsize=17)
            cbar2.ax.tick_params(labelsize=10)

            fig2.tight_layout()

            left = 0.1  # the left side of the subplots of the figure
            right = 0.90    # the right side of the subplots of the figure
            bottom = 0.09   # the bottom of the subplots of the figure
            top = 0.95     # the top of the subplots of the figure
            # wspace = 0.2  # the amount of width reserved for blank space between subplots
            # hspace = 0.2  # the amount of height reserved for white space between subplotss

            fig2.subplots_adjust(left=left, bottom=bottom,
                                 right=right, top=top)

            fig2.savefig(
                pwd + '/tg43appl_mbdca_compare_xz_isodose_profile_' + diameter + 'mmOut_'
                + vox_size + 'mm.pdf'
            )


def tg43_mbdca_comparison_histograms():
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots 
    from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' # for home
    target_dir = '/Users/student/research/results_not_to_git'  # for work

    air_kerma_true = 326.05715
    air_kerma_per_hist = 1.1584e-13
    max_dwell_time = 0.02917

    tg43_target_file = '/tg43appl_{1}mmOut_{0}pt0mm_sim.phantom_wo_applicator_wo_box.3ddose.gz'
    mlwa_target_file = '/mlwa_{1}mmOut_{0}shield_{2}pt0mm_sim.phantom_wo_applicator.3ddose.gz'
    base_title_str = 'TG43-MBDCA comparison \n {0}mm applicator; ncase = 5E9'

    outer_diams = [
        # '25',
        '30'
        # '35',
        # '40'
    ]

    shield_type_lst = [
        '0',
        '90',
        '180',
        '270',
    ]

    vox_size_lst = [
        '1',
        '2',
    ]

    dose_scale_factor = get_conversion_factor(
        air_kerma_true,
        air_kerma_per_hist,
        max_dwell_time
    )  # scale to maximum individual dwell time

    for diam_index, diameter in enumerate(outer_diams):

        for vox_index, vox_size in enumerate(vox_size_lst):

            title_str = base_title_str.format(diameter)

            fig, ax = subplots(
                2, 2, figsize=(10, 10),
                sharex='all', sharey='all'
            )

            tg43_full_data = DoseFile(
                target_dir + tg43_target_file.format(vox_size, diameter)
            )

            for shield_index, shield_type in enumerate(shield_type_lst):

                mlwa_full_data = DoseFile(
                    target_dir
                    + mlwa_target_file.format(shield_type, diameter, vox_size)
                )

                x_pos = array(tg43_full_data.positions[0])
                y_pos = array(tg43_full_data.positions[1])
                z_pos = array(tg43_full_data.positions[2])

                ax_x = shield_index // 2
                ax_y = shield_index % 2

                tg43_dose = tg43_full_data.dose * dose_scale_factor
                mlwa_dose = mlwa_full_data.dose * dose_scale_factor

                # calculate percentage difference between MBDCA calculation and tg43
                per_diff = nan_to_num(((tg43_dose - mlwa_dose) / mlwa_dose) * 100)

                xy_contour = ax[ax_x, ax_y].hist(
                    per_diff.flatten(),
                    bins='auto', density=True   
                )
                ax[ax_x, ax_y].set_title(
                    r'${0}$'.format(shield_type) + r'$^{\circ}$ shielding',
                    fontsize=20
                )

            for n in xrange(2):
                for m in xrange(2):

                    ax[n, m].grid(True)
                    ax[n, m].xaxis.set_tick_params(
                        labelsize=17, top=True, direction='in',rotation=45
                    )
                    ax[n, m].yaxis.set_tick_params(
                        labelsize=17, right=True, direction='in'
                    )
                    ax[n, m].set_xticks(arange(0, 1201, 150))
                    ax[n, m].set_xlim([-10, 1200])
                    ax[n, m].set_yticks(arange(0, 0.05 + 0.01, 0.01))
                    ax[n, m].set_ylim([0, 0.05])

            fig.text(
                0.02, 0.51, 'Probability density',
                fontsize=27, rotation='vertical', 
                va='center', ha='center'
            )
            fig.text(
                0.54, 0.025, 'Percent difference (%)', fontsize=27, 
                va='center', ha='center'
            )

            left = 0.12  # the left side of the subplots of the figure
            right = 0.96    # the right side of the subplots of the figure
            bottom = 0.12   # the bottom of the subplots of the figure
            top = 0.95     # the top of the subplots of the figure
            # wspace = 0.2  # the amount of width reserved for blank space between subplots
            # hspace = 0.2  # the amount of height reserved for white space between subplotss

            fig.tight_layout()
            fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

            fig.savefig(
                pwd + '/tg43-mbdca_compare_histograms_' + diameter + 'mmOut_'
                + vox_size + 'mm.pdf'
            )

if __name__ == "__main__":
    program_name, arguments = argv[0], argv[1:]

    # generate_tdvh_mlwa()
    # generate_tdvh_tg43()
    # dose_position_plots()
    # dose_inv_position_plots()
    # rel_dose_position_plot()
    # isodose_plot(mode='mlwa')
    isodose_plot_compare()
    # tg43_mbdca_comparison_isodose_plot()
    # tg43_mbdca_comparison_histograms()


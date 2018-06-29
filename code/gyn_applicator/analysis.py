#!/usr/bin/env
# author: Joseph Lucero
# created on: 3 May 2018 08:23:46
# purpose: plotting 3ddose files

from __future__ import division

from sys import argv, exit
from os import getcwd

from numpy import (
    linspace, zeros_like, histogram, arange, sqrt, isnan, array, meshgrid
    )
from scipy.interpolate import RegularGridInterpolator as RGI
from py3ddose import DoseFile, position_to_index
from normalize import get_conversion_factor

from matplotlib.cm import get_cmap
from matplotlib.style import use
use('seaborn-paper')
from matplotlib.pyplot import subplots

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

    # unshielded_file = target_dir + \
    #     '/mlwa_0shield_2pt0mm_sim.phantom_wo_applicator_wo_box.3ddose'

    # shielded_file_90 =  target_dir + \
    #     '/mlwa_90shield_2pt0mm_sim.phantom_wo_applicator_wo_box.3ddose'

    # shielded_file_180 = target_dir + \
    #     '/mlwa_180shield_2pt0mm_sim.phantom_wo_applicator_wo_box.3ddose'

    # shielded_file_270 = target_dir + \
    #     '/mlwa_270shield_2pt0mm_sim.phantom_wo_applicator_wo_box.3ddose'

    unshielded_file = target_dir + \
        '/mlwa_30mmOut_0shield_2pt0mm_sim.phantom_wo_applicator.3ddose.gz'

    shielded_file_90 =  target_dir + \
        '/mlwa_30mmOut_90shield_2pt0mm_sim.phantom_wo_applicator.3ddose.gz'

    shielded_file_180 = target_dir + \
        '/mlwa_30mmOut_180shield_2pt0mm_sim.phantom_wo_applicator.3ddose.gz'

    shielded_file_270 = target_dir + \
        '/mlwa_30mmOut_270shield_2pt0mm_sim.phantom_wo_applicator.3ddose.gz'

    title_str = 'Relative Dose vs. Position \n 192Ir-As-Source; ncase = 5E9'

    fig, ax = subplots(
        1, 1, figsize=(10, 10)
        )
    ax_twin = ax.twinx()

    unshielded_full_data = DoseFile(unshielded_file,load_uncertainty=True)
    shield_90_data = DoseFile(shielded_file_90,load_uncertainty=True)
    shield_180_data = DoseFile(shielded_file_180,load_uncertainty=True)
    shield_270_data = DoseFile(shielded_file_270,load_uncertainty=True)

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

    x_pos_mid = (x_pos[1:] + x_pos[:-1]) / 2.0
    y_pos_mid = (y_pos[1:] + y_pos[:-1]) / 2.0
    z_pos_mid = (z_pos[1:] + z_pos[:-1]) / 2.0

    x_position_to_index = {
        x_position: x_index
        for x_index, x_position in enumerate(unshielded_full_data.positions[0])
    }
    y_position_to_index = {
        y_position: y_index
        for y_index, y_position in enumerate(unshielded_full_data.positions[1])
    }
    z_position_to_index = {
        z_position: z_index
        for z_index, z_position in enumerate(unshielded_full_data.positions[2])
    }

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
        x_pos_mid[::2],
        # shield_90_data.dose[
        shield_90_interpolated_dose_matrix[
        # error_90[
            :, 
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_90[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$90^{\circ}$', color='darkgreen',
        markeredgecolor='darkgreen', marker='o',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )
    ax_twin.errorbar(
        x_pos_mid[::2],
        # shield_90_data.dose[
        shield_90_interpolated_dose_matrix[
        # error_90[
            :, 
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_90[
            :, 
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$90^{\circ}$', color='darkgreen',
        markeredgecolor='darkgreen', marker='o',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )

    ax.errorbar(
        x_pos_mid[::2],
        # shield_180_data.dose[
        shield_180_interpolated_dose_matrix[
        # error_180[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_180[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$180^{\circ}$', color='purple',
        markeredgecolor='purple', marker='s',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )
    ax_twin.errorbar(
        x_pos_mid[::2],
        # shield_180_data.dose[
        shield_180_interpolated_dose_matrix[
        # error_180[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_180[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$180^{\circ}$', color='purple',
        markeredgecolor='purple', marker='s',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )

    ax.errorbar(
        x_pos_mid[::2],
        # shield_270_data.dose[
        shield_270_interpolated_dose_matrix[
        # error_270[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_270[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$270^{\circ}$', color='darkorange',
        markeredgecolor='darkorange', marker='^',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )
    ax_twin.errorbar(
        x_pos_mid[::2],
        # shield_270_data.dose[
        shield_270_interpolated_dose_matrix[
        # error_270[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_270[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$270^{\circ}$', color='darkorange',
        markeredgecolor='darkorange', marker='^',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )

    ax.legend(loc=0, prop={'size': 18})
    ax.grid(False)
    ax.set_title('X - axis', fontsize=20)
    ax.set_yticks(arange(0.80, 1.11, 0.01))
    ax.set_ylim([0.80, 1.10])
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.set_xticks(arange(-10, 10 + 1, 1))
    ax.set_xlim([-10, 10])
    ax.axhline(1.0, 0, 0.5, lw=1.5, color='gray')

    ax_twin.grid(False)
    ax_twin.set_ylim([0.0, 0.3])
    ax_twin.set_yticks(arange(0.0, 0.31, 0.01))
    ax_twin.yaxis.set_tick_params(labelsize=14)
    ax_twin.vlines([-1.5,1.5],0,1.0)
    ax_twin.fill_between([-1.5, 1.5], 1.05, facecolor='lightgray')

    fig.text(
        0.02, 0.48, 'Relative Dose',
        fontsize=27, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.985, 0.51, 'Effective Transmission',
        fontsize=27, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.5, 0.51, 'A p p l i c a t o r',
        fontsize=35, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.43, 0.03, 'Position (cm)', fontsize=27, va='center'
    )
    fig.text(
        0.52, 0.95,
        title_str,
        fontsize=27, va='center', ha='center'
    )
    fig.tight_layout()

    left = 0.11  # the left side of the subplots of the figure
    right = 0.89    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    # wspace = 0.2  # the amount of width for blank space between subplots
    # hspace = 0.2  # the amount of height for white space between subplotss

    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig.savefig(pwd + '/relative_dosage_comparison.pdf')


def rel_dose_position_plot2():
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

    # unshielded_file = target_dir + \
    #     '/mlwa_0shield_2pt0mm_sim.phantom_wo_applicator_wo_box.3ddose'

    # shielded_file_90 =  target_dir + \
    #     '/mlwa_90shield_2pt0mm_sim.phantom_wo_applicator_wo_box.3ddose'

    # shielded_file_180 = target_dir + \
    #     '/mlwa_180shield_2pt0mm_sim.phantom_wo_applicator_wo_box.3ddose'

    # shielded_file_270 = target_dir + \
    #     '/mlwa_270shield_2pt0mm_sim.phantom_wo_applicator_wo_box.3ddose'

    unshielded_file = target_dir + \
        '/mlwa_30mmOut_0shield_2pt0mm_sim.phantom_wo_applicator.3ddose.gz'

    shielded_file_90 = target_dir + \
        '/mlwa_30mmOut_90shield_2pt0mm_sim.phantom_wo_applicator.3ddose.gz'

    shielded_file_180 = target_dir + \
        '/mlwa_30mmOut_180shield_2pt0mm_sim.phantom_wo_applicator.3ddose.gz'

    shielded_file_270 = target_dir + \
        '/mlwa_30mmOut_270shield_2pt0mm_sim.phantom_wo_applicator.3ddose.gz'

    title_str = 'Relative Dose vs. Position \n 192Ir-As-Source; ncase = 5E9'

    fig, ax = subplots(
        1, 1, figsize=(10, 10)
    )
    # ax_twin = ax.twinx()

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

    x_pos_mid = (x_pos[1:] + x_pos[:-1]) / 2.0
    y_pos_mid = (y_pos[1:] + y_pos[:-1]) / 2.0
    z_pos_mid = (z_pos[1:] + z_pos[:-1]) / 2.0

    x_position_to_index = {
        x_position: x_index
        for x_index, x_position in enumerate(unshielded_full_data.positions[0])
    }
    y_position_to_index = {
        y_position: y_index
        for y_index, y_position in enumerate(unshielded_full_data.positions[1])
    }
    z_position_to_index = {
        z_position: z_index
        for z_index, z_position in enumerate(unshielded_full_data.positions[2])
    }

    unshielded_interpolated_dose_matrix = unshielded_full_data.dose
    shield_90_interpolated_dose_matrix = shield_90_data.dose
    shield_180_interpolated_dose_matrix = shield_180_data.dose
    shield_270_interpolated_dose_matrix = shield_270_data.dose

    error_unshielded = unshielded_full_data.uncertainty * unshielded_interpolated_dose_matrix
    error_90 = shield_90_data.uncertainty * shield_90_interpolated_dose_matrix
    error_180 = shield_180_data.uncertainty * shield_180_interpolated_dose_matrix
    error_270 = shield_270_data.uncertainty * shield_270_interpolated_dose_matrix

    error_90 = sqrt(
        (error_90) ** 2
        + (error_unshielded) ** 2
    )
    error_180 = sqrt(
        (error_180) ** 2
        + (error_unshielded) ** 2
    )
    error_270 = sqrt(
        (error_270) ** 2
        + (error_unshielded) ** 2
    )

    error_90 = sqrt((error_90 / shield_90_interpolated_dose_matrix)**2 + (error_unshielded / unshielded_interpolated_dose_matrix)**2)
    error_180 = sqrt((error_180 / shield_180_interpolated_dose_matrix)**2 + (error_unshielded / unshielded_interpolated_dose_matrix)**2)
    error_270 = sqrt((error_270 / shield_270_interpolated_dose_matrix)**2 + (error_unshielded / unshielded_interpolated_dose_matrix)**2)

    shield_90_interpolated_dose_matrix -= unshielded_interpolated_dose_matrix
    shield_180_interpolated_dose_matrix -= unshielded_interpolated_dose_matrix
    shield_270_interpolated_dose_matrix -= unshielded_interpolated_dose_matrix

    shield_90_interpolated_dose_matrix /= unshielded_interpolated_dose_matrix
    shield_180_interpolated_dose_matrix /= unshielded_interpolated_dose_matrix
    shield_270_interpolated_dose_matrix /= unshielded_interpolated_dose_matrix

    shield_90_interpolated_dose_matrix = shield_90_interpolated_dose_matrix.__abs__() * 100
    shield_180_interpolated_dose_matrix = shield_180_interpolated_dose_matrix.__abs__() * 100
    shield_270_interpolated_dose_matrix = shield_270_interpolated_dose_matrix.__abs__() * 100

    error_90 *= shield_90_data.dose.__abs__()
    error_180 *= shield_180_data.dose.__abs__()
    error_270 *= shield_270_data.dose.__abs__()

    ax.errorbar(
        x_pos_mid[::2],
        # shield_90_data.dose[
        shield_90_interpolated_dose_matrix[
            # error_90[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_90[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$90^{\circ}$', color='darkgreen',
        markeredgecolor='darkgreen', marker='o',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )
    # ax_twin.errorbar(
    #     x_pos_mid[::2],
    #     # shield_90_data.dose[
    #     shield_90_interpolated_dose_matrix[
    #         # error_90[
    #         :,
    #         position_to_index(0.0, y_pos_mid),
    #         position_to_index(0.0, z_pos_mid)
    #     ][::2],
    #     yerr=error_90[
    #         :,
    #         position_to_index(0.0, y_pos_mid),
    #         position_to_index(0.0, z_pos_mid)
    #     ][::2],
    #     lw=0.0, label=r'$90^{\circ}$', color='darkgreen',
    #     markeredgecolor='darkgreen', marker='o',
    #     markeredgewidth=1, markersize=10, markerfacecolor='None',
    #     elinewidth=1.5, capsize=1.5
    # )

    ax.errorbar(
        x_pos_mid[::2],
        # shield_180_data.dose[
        shield_180_interpolated_dose_matrix[
            # error_180[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_180[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$180^{\circ}$', color='purple',
        markeredgecolor='purple', marker='s',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )
    # ax_twin.errorbar(
    #     x_pos_mid[::2],
    #     # shield_180_data.dose[
    #     shield_180_interpolated_dose_matrix[
    #         # error_180[
    #         :,
    #         position_to_index(0.0, y_pos_mid),
    #         position_to_index(0.0, z_pos_mid)
    #     ][::2],
    #     yerr=error_180[
    #         :,
    #         position_to_index(0.0, y_pos_mid),
    #         position_to_index(0.0, z_pos_mid)
    #     ][::2],
    #     lw=0.0, label=r'$180^{\circ}$', color='purple',
    #     markeredgecolor='purple', marker='s',
    #     markeredgewidth=1, markersize=10, markerfacecolor='None',
    #     elinewidth=1.5, capsize=1.5
    # )

    ax.errorbar(
        x_pos_mid[::2],
        # shield_270_data.dose[
        shield_270_interpolated_dose_matrix[
            # error_270[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        yerr=error_270[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::2],
        lw=0.0, label=r'$270^{\circ}$', color='darkorange',
        markeredgecolor='darkorange', marker='^',
        markeredgewidth=1, markersize=10, markerfacecolor='None',
        elinewidth=1.5, capsize=1.5
    )
    # ax_twin.errorbar(
    #     x_pos_mid[::2],
    #     # shield_270_data.dose[
    #     shield_270_interpolated_dose_matrix[
    #         # error_270[
    #         :,
    #         position_to_index(0.0, y_pos_mid),
    #         position_to_index(0.0, z_pos_mid)
    #     ][::2],
    #     yerr=error_270[
    #         :,
    #         position_to_index(0.0, y_pos_mid),
    #         position_to_index(0.0, z_pos_mid)
    #     ][::2],
    #     lw=0.0, label=r'$270^{\circ}$', color='darkorange',
    #     markeredgecolor='darkorange', marker='^',
    #     markeredgewidth=1, markersize=10, markerfacecolor='None',
    #     elinewidth=1.5, capsize=1.5
    # )

    ax.legend(loc=0, prop={'size': 18})
    ax.grid(True)
    ax.set_title('X - axis', fontsize=20)
    # ax.set_yticks(arange(0.80, 1.11, 0.01))
    ax.set_ylim([0.0, 100.0])
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.set_xticks(arange(-10, 10 + 1, 1))
    ax.set_xlim([-10, 10])
    # ax.axhline(1.0, 0, 0.5, lw=1.5, color='gray')

    # ax_twin.grid(False)
    # ax_twin.set_ylim([0.0, 0.3])
    # ax_twin.set_yticks(arange(0.0, 0.31, 0.01))
    # ax_twin.yaxis.set_tick_params(labelsize=14)
    # ax_twin.vlines([-1.5, 1.5], 0, 1.0)
    # ax_twin.fill_between([-1.5, 1.5], 1.05, facecolor='lightgray')

    fig.text(
        0.02, 0.48, 'Relative Dose',
        fontsize=27, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.985, 0.51, 'Effective Transmission',
        fontsize=27, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.5, 0.51, 'A p p l i c a t o r',
        fontsize=35, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.43, 0.03, 'Position (cm)', fontsize=27, va='center'
    )
    fig.text(
        0.52, 0.95,
        title_str,
        fontsize=27, va='center', ha='center'
    )
    fig.tight_layout()

    left = 0.11  # the left side of the subplots of the figure
    right = 0.89    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    # wspace = 0.2  # the amount of width for blank space between subplots
    # hspace = 0.2  # the amount of height for white space between subplotss

    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig.savefig(pwd + '/relative_dosage_comparison.pdf')

def isodose_plot():
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

    file_list = [
        target_dir +
        '/mlwa_30mmOut_0shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        target_dir +
        '/mlwa_30mmOut_90shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        target_dir +
        '/mlwa_30mmOut_180shield_{0}_sim.phantom_wo_applicator.3ddose.gz',
        target_dir +
        '/mlwa_30mmOut_270shield_{0}_sim.phantom_wo_applicator.3ddose.gz'
        # target_dir +
        # '/mlwa_0shield_{0}_sim.phantom_wo_applicator.3ddose',
        # target_dir +
        # '/mlwa_90shield_{0}_sim.phantom_wo_applicator.3ddose',
        # target_dir +
        # '/mlwa_180shield_{0}_sim.phantom_wo_applicator.3ddose',
        # target_dir +
        # '/mlwa_270shield_{0}_sim.phantom_wo_applicator.3ddose'
    ]

    shield_type_lst = [
        'Unshielded',
        r'$90^{\circ}$ shield',
        r'$180^{\circ}$ shield',
        r'$270^{\circ}$ shield'
    ]

    vox_size_lst = [
        # '1mm',
        '2mm'
    ]

    vox_size_txt_lst = [
        # '1pt0mm',
        '2pt0mm'
    ]

    # title_str = 'Isodose Contours \n Applicator-As-Source; ncase = 5E9'
    title_str = 'Isodose Contours \n 192Ir-As-Source; ncase = 5E9'

    for fig_index in xrange(len(vox_size_lst)):

        fig, ax = subplots(
            2, 2, figsize=(12, 12),
            sharex='all', sharey='all'
        )
        fig2, ax2 = subplots(
            2, 2, figsize=(12, 12),
            sharex='all', sharey='all'
        )

        for file_index, file in enumerate(file_list):

            full_data = DoseFile(file.format(vox_size_txt_lst[fig_index]))

            y_position_to_index = {
                position: index
                for index, position in enumerate(full_data.positions[1])
            }
            z_position_to_index = {
                position: index
                for index, position in enumerate(full_data.positions[2])
            }

            x_pos = array(full_data.positions[0])
            y_pos = array(full_data.positions[1])
            z_pos = array(full_data.positions[2])

            x_pos_mid = ((x_pos[1:]) + (x_pos[:-1])) / 2.0
            y_pos_mid = ((y_pos[1:]) + (y_pos[:-1])) / 2.0
            z_pos_mid = ((z_pos[1:]) + (z_pos[:-1])) / 2.0

            if file_index == 0:
                ax_x, ax_y = 0, 0
            elif file_index == 1:
                ax_x, ax_y = 0, 1
            elif file_index == 2:
                ax_x, ax_y = 1, 0
            else:
                ax_x, ax_y = 1, 1

            Nx, Ny, Nz = full_data.shape

            full_data.dose *= 8.2573429808917e13  # scale to maximum individual dwell time

            full_data.dose /= 5  # normalize to desired dose of 5 Gy
            full_data.dose *= 100  # express in percent. Should see 100% at x=-2cm

            x_min, x_max = full_data.x_extent
            y_min, y_max = full_data.y_extent
            z_min, z_max = full_data.z_extent

            xy_contour = ax[ax_x, ax_y].contourf(
                x_pos_mid, y_pos_mid,
                # matplotlib plots column by row (instead of row by column)
                # so transpose data array to account for this
                full_data.dose[:, :, position_to_index(
                    0.0, z_pos_mid)].transpose(),
                arange(0, 110, 10),
                # [5, 10, 20, 50, 100],
                # cmap=get_cmap('gnuplot')
            )
            # ax[ax_x, ax_y].contour(
            #     x_pos_mid, y_pos_mid,
            #     # matplotlib plots column by row (instead of row by column)
            #     # so transpose data array to account for this
            #     full_data.dose[:, :, position_to_index(0.0, z_pos_mid)].transpose(),
            #     # arange(0, 110, 10),
            #     [5, 10, 20, 50, 100],
            #     cmap=get_cmap('RdPu')
            # )
            ax[ax_x, ax_y].set_title(shield_type_lst[file_index], fontsize=15)

            xz_contour = ax2[ax_x, ax_y].contourf(
                x_pos_mid, z_pos_mid,
                # matplotlib plots column by row (instead of row by column)
                # so transpose data array to account for this
                full_data.dose[:, position_to_index(
                    0.0, y_pos_mid), :].transpose(),
                arange(0, 110, 10),
                # [5, 10, 20, 50, 100],
                # cmap=get_cmap('gnuplot')
            )
            # ax2[ax_x, ax_y].contour(
            #     x_pos_mid, z_pos_mid,
            #     # matplotlib plots column by row (instead of row by column)
            #     # so transpose data array to account for this
            #     full_data.dose[:, position_to_index(0.0, y_pos_mid), :].transpose(),
            #     # arange(0, 110, 10),
            #     [5, 10, 20, 50, 100],
            #     cmap=get_cmap('RdPu')
            # )
            ax2[ax_x, ax_y].set_title(shield_type_lst[file_index], fontsize=15)

        for n in xrange(2):
            for m in xrange(2):

                # ax[n, m].grid(True)
                ax[n, m].xaxis.set_tick_params(
                    labelsize=14, top=True, direction='in'
                )
                ax[n, m].yaxis.set_tick_params(
                    labelsize=14, right=True, direction='in'
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
                    labelsize=14, top=True, direction='in'
                )
                ax2[n, m].yaxis.set_tick_params(
                    labelsize=14, right=True, direction='in'
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
        fig.text(
            0.52, 0.95,
            title_str,
            fontsize=27, va='center', ha='center'
        )

        cax = fig.add_axes([0.91, 0.13, 0.01, 0.7])
        cbar1 = fig.colorbar(
            xy_contour, cax=cax, orientation='vertical',
            ax=ax
        )
        cbar1.set_label('Percentage Isodose (%)', fontsize=24)
        cbar1.ax.tick_params(labelsize=14)
        fig.tight_layout()

        left = 0.1  # the left side of the subplots of the figure
        right = 0.89    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.88     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width reserved for blank space between subplots
        # hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig.savefig(
            pwd + '/xy_isodose_profile_' + vox_size_lst[fig_index] + '.pdf'
        )

        fig2.text(
            0.01, 0.51, 'z-axis (cm)',
            fontsize=27, rotation='vertical', va='center'
        )
        fig2.text(
            0.43, 0.03, 'x-axis (cm)', fontsize=27, va='center'
        )
        fig2.text(
            0.52, 0.95,
            title_str,
            fontsize=27, va='center', ha='center'
        )
        cax2 = fig2.add_axes([0.91, 0.13, 0.01, 0.7])
        cbar2 = fig2.colorbar(
            xz_contour, cax=cax2, orientation='vertical',
            ax=ax
        )
        cbar2.set_label('Percentage Isodose (%)', fontsize=24)
        cbar2.ax.tick_params(labelsize=14)

        fig2.tight_layout()

        left = 0.1  # the left side of the subplots of the figure
        right = 0.89    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.88     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width reserved for blank space between subplots
        # hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig2.savefig(
            pwd + '/xz_isodose_profile_' + vox_size_lst[fig_index] + '.pdf'
        )
    
if __name__ == "__main__":
    program_name, arguments = argv[0], argv[1:]

    # generate_tdvh_mlwa()
    # generate_tdvh_tg43()
    # dose_position_plots()
    # rel_dose_position_plot()
    rel_dose_position_plot2()
    # isodose_plot()


#!/usr/bin/env
# author: Joseph Lucero
# created on: 3 May 2018 08:23:46
# purpose: plotting 3ddose files

from __future__ import division

from sys import argv, exit
from os import getcwd

from numpy import linspace, zeros_like, histogram, arange 
from py3ddose import DoseFile

from matplotlib import cm
from matplotlib.style import use
use('seaborn')
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
            # data_base_flat = main_base_data.dose.flatten() * 2.2861e14 # converts to Gy; norm to treatment time
            data_flat = main_data.dose.flatten() * 8.2573429808917e13 # converts to Gy; norm to individual max dwell time
            
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
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

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
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig2.subplots_adjust(left=left, top=top, right=right, bottom=bottom)
    fig2.savefig(pwd + '/cumulative_dose_volume_histogram.pdf')


def generate_tdvh_tg43():

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' # for home
    target_dir = '/Users/student/research/results_not_to_git'  # for work

    file_tg43_0point5mm = target_dir + \
        '/tg43_0shield_0.5mm_sim.phantom_wo_box.3ddose'
    file_tg43_1point0mm = target_dir + \
        '/tg43_0shield_1mm_sim.phantom_wo_box.3ddose'
    file_tg43_2point0mm = target_dir + \
        '/tg43_0shield_2mm_sim.phantom_wo_box.3ddose'

    file_tg43_applicator_0point5mm = target_dir + \
        '/tg43_applicator_0.5mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_tg43_applicator_1point0mm = target_dir + \
        '/tg43_applicator_1mm_sim.phantom_wo_applicator_wo_box.3ddose'
    file_tg43_applicator_2point0mm = target_dir + \
        '/tg43_applicator_2mm_sim.phantom_wo_applicator_wo_box.3ddose'

    file_dict = {
        (0, 0): file_tg43_0point5mm,
        (0, 1): file_tg43_1point0mm,
        (0, 2): file_tg43_2point0mm,

        (1, 0): file_tg43_applicator_0point5mm,
        (1, 1): file_tg43_applicator_1point0mm,
        (1, 2): file_tg43_applicator_2point0mm
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

    fig, ax = subplots(2, 1, figsize=(10, 10), sharex='all', sharey='all')
    fig2, ax2 = subplots(2, 1, figsize=(10, 10), sharex='all', sharey='all')

    for index1 in xrange(2):
        for index2 in xrange(3):

            main_data = DoseFile(file_dict[(index1, index2)])
            # data_base_flat = main_base_data.dose.flatten() * 2.2861e14 # converts to Gy; norm to treatment time
            # converts to Gy; norm to individual max dwell time
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
        ax[i].legend(loc=0, prop={'size': 10})
        ax[i].set_title(title_list[i], fontsize=20)
        ax[i].xaxis.set_tick_params(labelsize=17)
        ax[i].yaxis.set_tick_params(labelsize=17)

        ax2[i].legend(loc=0, prop={'size': 10})
        ax2[i].set_title(title_list[i], fontsize=20)
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
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

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
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

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
        '/mlwa_0shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_90shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_180shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_270shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
    ]

    label_list = [
        'Unshielded',
        '90 degree shield',
        '180 degree shield',
        '270 degree shield'
    ]

    vox_size_list = [
        '0.5mm',
        '1mm',
        '2mm',
    ]

    color_list = [
        'darkgreen',
        'purple',
        'darkorange',
        'saddlebrown'
    ]

    title_list = [
        'X - axis',
        'Y - axis',
        'Z - axis'
    ]

    for fig_index in xrange(3):

        fig, ax = subplots(
            3, 1, figsize=(10, 10),
            sharex='col'
        )

        for index2 in xrange(4):  # iterate through shield types

            full_data = DoseFile(
                file_list[index2].format(vox_size_list[fig_index])
            )

            Nx, Ny, Nz = full_data.shape

            # full_data.dose *= 2.2861e14 # scale to total treatment time
            full_data.dose *= 8.2573429808917e13  # scale to maximum individual dwell time

            x_min, x_max = full_data.x_extent
            y_min, y_max = full_data.y_extent
            z_min, z_max = full_data.z_extent

            for index1 in xrange(3):  # iterate through axes

                if index1 == 0:
                    ax[index1, index2].plot(
                        linspace(x_min, x_max, Nx), full_data.dose[Nz // 2, Ny // 2, :],
                        # yerr=full_data.uncertainty[Nz // 2, Ny // 2, :],
                        lw=3.0
                    )
                elif index1 == 1:
                    ax[index1, index2].plot(
                        linspace(y_min, y_max, Ny), full_data.dose[Nz // 2, :, Nx // 2],
                        # yerr=full_data.uncertainty[Nz // 2, :, Nx // 2],
                        lw=3.0
                    )
                else:
                    ax[index1, index2].plot(
                        linspace(z_min, z_max, Nz), full_data.dose[:, Ny // 2, Nx // 2],
                        # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                        lw=3.0
                    )

        for n in xrange(3):
            for m in xrange(4):
                # ax[n].legend(loc=0,prop={'size':8})
                ax[n, m].grid(True)
                if n == 0:
                    ax[n, m].set_title(title_list[n],fontsize=20)
                ax[n, m].xaxis.set_tick_params(labelsize=14)
                ax[n, m].yaxis.set_tick_params(labelsize=14)
                if n == 2:
                    ax[n, m].set_xticks(arange(z_min, z_max))

        fig.text(
            0.01, 0.51, 'Dose (Gy)',
            fontsize=27, rotation='vertical', va='center'
        )
        fig.text(
            0.43, 0.03, 'Position (cm)', fontsize=27, va='center'
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
        wspace = 0.2  # the amount of width reserved for blank space between subplots
        hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig.savefig(pwd + '/dosage_comparison_' + vox_size_list[fig_index] + '.pdf')


def rel_dose_position_plot():
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

    unshielded_file = target_dir + \
        '/mlwa.phantom_wo_applicator.3ddose'

    shielded_file_90 =  target_dir + \
        '/mlwa_90.phantom_wo_applicator.3ddose'

    shielded_file_180 = target_dir + \
        '/mlwa_180.phantom_wo_applicator.3ddose'

    shielded_file_270 = target_dir + \
        '/mlwa_270.phantom_wo_applicator.3ddose',

    fig, ax = subplots(
        1, 1, figsize=(10, 10),
        sharex='col'
    )
    ax_twin = ax.twinx()

    unshielded_full_data = DoseFile(unshielded_file)
    shield_90_data = DoseFile(shielded_file_90)
    shield_180_data = DoseFile(shielded_file_180)
    shield_270_data = DoseFile(shielded_file_270)

    Nx, Ny, Nz = unshielded_full_data.shape

    shield_90_data.dose /= unshielded_full_data.dose
    shield_180_data.dose /= unshielded_full_data.dose
    shield_270_data.dose /= unshielded_full_data.dose

    x_min, x_max = unshielded_full_data.x_extent
    y_min, y_max = unshielded_full_data.y_extent
    z_min, z_max = unshielded_full_data.z_extent

    ax.plot(
        linspace(x_min, -1.6, (Nx//2 - 8) // 2),
        shield_90_data.dose[Nz // 2, Ny // 2, :Nx//2 - 8][::2],
        lw=0.0, label='90 degrees',
        color='darkgreen', marker='o',
        markersize=10
    )
    ax_twin.plot(
        linspace(1.6, x_max, (Nx//2 - 8) // 2),
        shield_90_data.dose[Nz // 2, Ny // 2, Nx//2 + 8:][::2],
        lw=0.0, color='darkgreen', marker='o',
        markersize=10
    )

    ax.plot(
        linspace(x_min, -1.6, (Nx//2 - 8) // 2),
        shield_180_data.dose[Nz // 2, Ny // 2, :Nx//2 - 8][::2],
        lw=0.0, label='180 degrees',
        color='purple', marker='s',
        markersize=10
    )
    ax_twin.plot(
        linspace(1.6, x_max, (Nx//2 - 8) // 2),
        shield_180_data.dose[Nz // 2, Ny // 2, Nx//2 + 8:][::2],
        lw=0.0, color='purple', marker='s',
        markersize=10
    )

    ax.plot(
        linspace(x_min, -1.6, (Nx//2 - 8) // 2),
        shield_270_data.dose[Nz // 2, Ny // 2, :Nx//2 - 8][::2],
        lw=0.0, label='270 degrees',
        color='darkorange', marker='^',
        markersize=10
    )
    ax_twin.plot(
        linspace(1.6, x_max, (Nx//2 - 8) // 2),
        shield_270_data.dose[Nz // 2, Ny // 2, Nx//2 + 8:][::2],
        lw=0.0, color='darkorange', marker='^',
        markersize=10
    )

    ax.legend(loc=0, prop={'size': 18})
    ax.grid(True)
    ax.set_title('X - axis', fontsize=20)
    ax.set_ylim([0.8, 1.05])
    ax.xaxis.set_tick_params(labelsize=17)
    ax.yaxis.set_tick_params(labelsize=17)
    ax_twin.set_ylim([0.0, 0.5])
    ax_twin.yaxis.set_tick_params(labelsize=17)
    ax_twin.vlines([-1.5,1.5],0,1.0)

    ax.set_xticks(arange(z_min, z_max + 1, 2))

    fig.text(
        0.01, 0.48, 'Relative Dose',
        fontsize=27, rotation='vertical', va='center'
    )
    fig.text(
        0.96, 0.51, 'Effective Transmission',
        fontsize=27, rotation='vertical', va='center'
    )
    fig.text(
        0.43, 0.03, 'Position (cm)', fontsize=27, va='center'
    )
    fig.text(
        0.52, 0.95,
        'Relative Dose vs. Position \n With Volume Correction; ncase = 1E9',
        fontsize=27, va='center', ha='center'
    )
    fig.tight_layout()

    left = 0.11  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

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
        '/mlwa_0shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_90shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_180shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_270shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
    ]

    label_list = [
        'Unshielded',
        '90 degree shield',
        '180 degree shield',
        '270 degree shield'
    ]

    vox_size_list = [
        '0.5mm',
        '1mm',
        '2mm',
    ]

    for fig_index in xrange(3):

        fig, ax = subplots(
            2, 2, figsize=(10, 10),
            sharex='all', sharey='all'
        )
        fig2, ax2 = subplots(
            2, 2, figsize=(10, 10),
            sharex='all', sharey='all'
        )

        for file_index, file in enumerate(file_list):

            full_data = DoseFile(file.format(vox_size_list[fig_index]))

            if file_index == 0:
                ax_x, ax_y = 0, 0
            elif file_index == 1:
                ax_x, ax_y = 0, 1
            elif file_index == 2:
                ax_x, ax_y = 1, 0
            else:
                ax_x, ax_y = 1, 1

            Nx, Ny, Nz = full_data.shape

            # full_data.dose *= 2.2861e14 # scale to total treatment time
            full_data.dose *= 8.2573429808917e13  # scale to maximum individual dwell time

            full_data.dose /= 5  # normalize to desired dose of 5 Gy
            full_data.dose *= 100  # express in percent. Should see 100% at x=-2cm

            x_min, x_max = full_data.x_extent
            y_min, y_max = full_data.y_extent
            z_min, z_max = full_data.z_extent

            xy_contour = ax[ax_x, ax_y].contourf(
                linspace(x_min, x_max, Nx), linspace(y_min, y_max, Ny),
                full_data.dose[Nz // 2, :, :],
                arange(0, 110, 10),
                cmap=cm.Purples
            )

            xz_contour = ax2[ax_x, ax_y].contourf(
                linspace(x_min, x_max, Nx), linspace(z_min, z_max, Nz),
                full_data.dose[:, Ny // 2, :],
                arange(0, 110, 10),
                cmap=cm.Purples
            )

        for n in xrange(2): 
            for m in xrange(2):

                ax[n, m].grid(True)
                ax[n, m].xaxis.set_tick_params(labelsize=14)
                ax[n, m].yaxis.set_tick_params(labelsize=14)
                ax[n, m].set_xticks(arange(x_min, x_max + 1, 2))
                ax[n, m].set_yticks(arange(y_min, y_max + 1, 2))

                ax2[n ,m].xaxis.set_tick_params(labelsize=14)
                ax2[n ,m].yaxis.set_tick_params(labelsize=14)
                ax2[n ,m].set_xticks(arange(x_min, x_max + 1, 2))
                ax2[n ,m].set_yticks(arange(y_min, y_max + 1, 2))
                ax2[n ,m].vlines([-2, 2],-10, 10,linestyles='dashed',lw=2.0)
                ax2[n ,m].grid(True)

        fig.text(
            0.01, 0.51, 'y-axis (cm)',
            fontsize=27, rotation='vertical', va='center'
        )
        fig.text(
            0.43, 0.03, 'x-axis (cm)', fontsize=27, va='center'
        )
        fig.text(
            0.52, 0.95,
            'Isodose Contours \n With Volume Correction; ncase = 1E9',
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
        top = 0.90     # the top of the subplots of the figure
        wspace = 0.2  # the amount of width reserved for blank space between subplots
        hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig.savefig(
            pwd + '/xy_isodose_profile_' + vox_size_list[fig_index] + '.pdf'
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
            'Isodose Contours \n With Volume Correction; ncase = 1E9',
            fontsize=27, va='center', ha='center'
        )
        cax2 = fig2.add_axes([0.91, 0.13, 0.01, 0.7])
        cbar2 = fig2.colorbar(
            xz_contour, cax=cax2, orientation='vertical',
            ax=ax
        )
        cbar2.set_label('Percentage Isodose (%)', fontsize=24)
        cbar2.set_clim([0, 100])
        cbar2.ax.tick_params(labelsize=14)

        fig2.tight_layout()

        left = 0.1  # the left side of the subplots of the figure
        right = 0.89    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.90     # the top of the subplots of the figure
        wspace = 0.2  # the amount of width reserved for blank space between subplots
        hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig2.savefig(
            pwd + '/xz_isodose_profile_' + vox_size_list[fig_index] + '.pdf'
            )
    
if __name__ == "__main__":
    program_name, arguments = argv[0], argv[1:]

    # generate_tdvh_mlwa()
    # generate_tdvh_tg43()
    # dose_position_plots()
    # rel_dose_position_plot()
    isodose_plot()


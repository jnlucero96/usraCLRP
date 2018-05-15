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

def generate_tdvh():

    pwd = getcwd()

    # target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' # for home
    target_dir = '/Users/student/research/results_not_to_git' # for work

    file_base_dict = {
        0:target_dir + '/mlwa.phantom_wo_applicator.3ddose',
        1:target_dir + '/mlwa_2.phantom_wo_applicator.3ddose',
        2:target_dir + '/mlwa_3.phantom_wo_applicator.3ddose',
        3:target_dir + '/mlwa_4.phantom_wo_applicator.3ddose',
        4:target_dir + '/mlwa_5.phantom_wo_applicator.3ddose',
        5:target_dir + '/mlwa_6.phantom_wo_applicator.3ddose',
        6:target_dir + '/mlwa_7.phantom_wo_applicator.3ddose',
        7:target_dir + '/mlwa_8.phantom_wo_applicator.3ddose'
    }
    
    file_90_dict = {
        0:target_dir + '/mlwa_90.phantom_wo_applicator.3ddose',
        1:target_dir + '/mlwa_90_2.phantom_wo_applicator.3ddose',
        2:target_dir + '/mlwa_90_3.phantom_wo_applicator.3ddose',
        3:target_dir + '/mlwa_90_4.phantom_wo_applicator.3ddose',
        4:target_dir + '/mlwa_90_5.phantom_wo_applicator.3ddose',
        5:target_dir + '/mlwa_90_6.phantom_wo_applicator.3ddose',
        6:target_dir + '/mlwa_90_7.phantom_wo_applicator.3ddose',
        7:target_dir + '/mlwa_90_8.phantom_wo_applicator.3ddose'
    }
    
    file_180_dict = {
        0:target_dir + '/mlwa_180.phantom_wo_applicator.3ddose',
        1:target_dir + '/mlwa_180_2.phantom_wo_applicator.3ddose',
        2:target_dir + '/mlwa_180_3.phantom_wo_applicator.3ddose',
        3:target_dir + '/mlwa_180_4.phantom_wo_applicator.3ddose',
        4:target_dir + '/mlwa_180_5.phantom_wo_applicator.3ddose',
        5:target_dir + '/mlwa_180_6.phantom_wo_applicator.3ddose',
        6:target_dir + '/mlwa_180_7.phantom_wo_applicator.3ddose',
        7:target_dir + '/mlwa_180_8.phantom_wo_applicator.3ddose'
    }
    
    file_270_dict = {
        0:target_dir + '/mlwa_270.phantom_wo_applicator.3ddose',
        1:target_dir + '/mlwa_270_2.phantom_wo_applicator.3ddose',
        2:target_dir + '/mlwa_270_3.phantom_wo_applicator.3ddose',
        3:target_dir + '/mlwa_270_4.phantom_wo_applicator.3ddose',
        4:target_dir + '/mlwa_270_5.phantom_wo_applicator.3ddose',
        5:target_dir + '/mlwa_270_6.phantom_wo_applicator.3ddose',
        6:target_dir + '/mlwa_270_7.phantom_wo_applicator.3ddose',
        7:target_dir + '/mlwa_270_8.phantom_wo_applicator.3ddose'
    }

    label_list = [
        'Base Config.',
        'Config. 2',
        'Config. 3',
        'Config. 4',
        'Config. 5',
        'Config. 6',
        'Config. 7',
        'Config. 8'
    ]

    color_list = [
        'darkgreen',
        'purple',
        'darkorange',
        'red',
        'black',
        'blue',
        'saddlebrown',
        'darkcyan'
    ]

    title_list = [
        'Unshielded',
        '90 Degree Shielding',
        '180 Degree Shielding',
        '270 Degree Shielding'
    ]

    fig, ax = subplots(4, 1, figsize=(10, 10), sharex='all', sharey='all')
    fig2, ax2 = subplots(4, 1, figsize=(10, 10), sharex='all', sharey='all')

    for index in xrange(8):

        main_base_data = DoseFile(file_base_dict[index])
        # data_base_flat = main_base_data.dose.flatten() * 2.2861e14 # converts to Gy; norm to treatment time
        data_base_flat = main_base_data.dose.flatten() * 8.2573429808917e13 # converts to Gy; norm to individual max dwell time
        
        main_90_data = DoseFile(file_90_dict[index])
        # data_90_flat = main_90_data.dose.flatten() * 2.2861e14 # converts to Gy; norm to treatment time
        data_90_flat = main_90_data.dose.flatten() * 8.2573429808917e13 # converts to Gy; norm to individual max dwell time
        
        main_180_data = DoseFile(file_180_dict[index])
        # data_180_flat = main_180_data.dose.flatten() * 2.2861e14 # converts to Gy; norm to treatment time
        data_180_flat = main_180_data.dose.flatten() * 8.2573429808917e13 # converts to Gy; norm to individual max dwell time

        main_270_data = DoseFile(file_270_dict[index])
        # data_270_flat = main_270_data.dose.flatten() * 2.2861e14  # converts to Gy; norm to treatment time
        data_270_flat = main_270_data.dose.flatten() * 8.2573429808917e13  # converts to Gy; norm to individual max dwell time

        n_base, bin_base, __ = ax[0].hist(
            data_base_flat,
            bins=500,
            color=color_list[index],
            label=label_list[index],
            weights=zeros_like(data_base_flat) + 1. / data_base_flat.size * 100,
            alpha=0.4
            )
        
        n_90, bin_90, __ = ax[1].hist(
            data_90_flat,
            bins=500,
            color=color_list[index],
            label=label_list[index],
            weights=zeros_like(data_90_flat) + 1. / data_90_flat.size * 100,
            alpha=0.4
            )
        
        n_180, bin_180, __ = ax[2].hist(
            data_180_flat,
            bins=500,
            color=color_list[index],
            label=label_list[index],
            weights=zeros_like(data_180_flat) + 1. / data_180_flat.size * 100,
            alpha=0.4
            )
            
        n_270, bin_270, __ = ax[3].hist(
            data_270_flat,
            bins=500,
            color=color_list[index],
            label=label_list[index],
            weights=zeros_like(data_270_flat) + 1. / data_270_flat.size * 100,
            alpha=0.4
            )

        n_cum_base = n_base[::-1].cumsum()[::-1]
        n_cum_90 = n_90[::-1].cumsum()[::-1]
        n_cum_180 = n_180[::-1].cumsum()[::-1]
        n_cum_270 = n_270[::-1].cumsum()[::-1]

        ax2[0].loglog(
            bin_base[:-1], n_cum_base,
            color=color_list[index],
            label=label_list[index]
        )
        ax2[1].loglog(
            bin_90[:-1], n_cum_90,
            color=color_list[index],
            label=label_list[index]
        )
        ax2[2].loglog(
            bin_180[:-1], n_cum_180,
            color=color_list[index],
            label=label_list[index]
        )
        ax2[3].loglog(
            bin_270[:-1], n_cum_270,
            color=color_list[index],
            label=label_list[index]
        )

    for i in xrange(4):
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
        target_dir + '/mlwa.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_2.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_3.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_4.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_5.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_6.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_7.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_8.phantom_wo_applicator.3ddose'
    ]

    label_list = [
        'Base Config.',
        'Config. 2',
        'Config. 3',
        'Config. 4',
        'Config. 5',
        'Config. 6',
        'Config. 7',
        'Config. 8'
    ]

    color_list = [
        'darkgreen',
        'purple',
        'darkorange',
        'red',
        'black',
        'blue',
        'saddlebrown',
        'darkcyan'
    ]

    title_list = [
        'X - axis',
        'Y - axis',
        'Z - axis'
    ]

    fig, ax = subplots(
        3, 1, figsize=(10, 10),
        sharex='col'
    )

    for index, file_name in enumerate(file_list):

        full_data = DoseFile(file_name)

        Nx, Ny, Nz = full_data.shape

        # full_data.dose *= 2.2861e14 # scale to total treatment time
        full_data.dose *= 8.2573429808917e13  # scale to maximum individual dwell time

        x_min, x_max = full_data.x_extent
        y_min, y_max = full_data.y_extent
        z_min, z_max = full_data.z_extent

        ax[0].plot(
            linspace(x_min, x_max, Nx), full_data.dose[Nz // 2, Ny // 2, :],
            # yerr=full_data.uncertainty[Nz // 2, Ny // 2, :],
            lw=3.0, label=label_list[index], color=color_list[index]
        )
        ax[1].plot(
            linspace(y_min, y_max, Ny), full_data.dose[Nz // 2, :, Nx // 2],
            # yerr=full_data.uncertainty[Nz // 2, :, Nx // 2],
            lw=3.0, label=label_list[index], color=color_list[index]
        )
        ax[2].plot(
            linspace(z_min, z_max, Nz), full_data.dose[:, Ny // 2, Nx // 2],
            # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
            lw=3.0, label=label_list[index], color=color_list[index]
        )

    for n in xrange(3):
        ax[n].legend(loc=0,prop={'size':8})
        ax[n].grid(True)
        ax[n].set_title(title_list[n],fontsize=20)
        ax[n].xaxis.set_tick_params(labelsize=17)
        ax[n].yaxis.set_tick_params(labelsize=17)

    ax[2].set_xticks(arange(z_min,z_max))

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

    left = 0.1  # the left side of the subplots of the figure
    right = 0.95    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig.savefig(pwd + '/dosage_comparison.pdf')


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

    unshielded_file_list = [
        target_dir + '/mlwa.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_2.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_3.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_4.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_5.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_6.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_7.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_8.phantom_wo_applicator.3ddose'
    ]

    shielded_file_90_list = [
        target_dir + '/mlwa_90.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_90_2.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_90_3.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_90_4.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_90_5.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_90_6.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_90_7.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_90_8.phantom_wo_applicator.3ddose'
    ]

    shielded_file_180_list = [
        target_dir + '/mlwa_180.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_180_2.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_180_3.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_180_4.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_180_5.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_180_6.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_180_7.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_180_8.phantom_wo_applicator.3ddose'
    ]

    shielded_file_270_list = [
        target_dir + '/mlwa_270.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_270_2.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_270_3.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_270_4.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_270_5.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_270_6.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_270_7.phantom_wo_applicator.3ddose',
        target_dir + '/mlwa_270_8.phantom_wo_applicator.3ddose'
    ]

    label_list = [
        'Base Config.',
        'Config. 2',
        'Config. 3',
        'Config. 4',
        'Config. 5',
        'Config. 6',
        'Config. 7',
        'Config. 8'
    ]

    color_list = [
        'darkgreen',
        'purple',
        'darkorange',
        'red',
        'black',
        'blue',
        'saddlebrown',
        'darkcyan'
    ]

    title_list = [
        'X - axis',
        'Y - axis',
        'Z - axis'
    ]

    fig, ax = subplots(
        1, 1, figsize=(10, 10),
        sharex='col'
    )
    ax_twin = ax.twinx()

    for index, file_90_name in enumerate(shielded_file_90_list):

        if index == 1:
            

            unshielded_full_data = DoseFile(unshielded_file_list[index])
            shield_90_data = DoseFile(file_90_name)
            shield_180_data = DoseFile(shielded_file_180_list[index])
            shield_270_data = DoseFile(shielded_file_270_list[index])

            Nx, Ny, Nz = unshielded_full_data.shape

            # full_data.dose *= 2.2861e14 # scale to total treatment time
            # full_data.dose *= 8.2573429808917e13  # scale to maximum individual dwell time

            shield_90_data.dose /= unshielded_full_data.dose
            shield_180_data.dose /= unshielded_full_data.dose
            shield_270_data.dose /= unshielded_full_data.dose

            x_min, x_max = unshielded_full_data.x_extent
            y_min, y_max = unshielded_full_data.y_extent
            z_min, z_max = unshielded_full_data.z_extent

            ax.plot(
                linspace(x_min, -1.6, (Nx//2 - 8) // 2), 
                shield_90_data.dose[Nz // 2 + 25, Ny // 2, :Nx//2 - 8][::2],
                lw=0.0, label='90 degrees',
                color='darkgreen', marker='o', 
                markersize=10
            )
            ax_twin.plot(
                linspace(1.6, x_max, (Nx//2 - 8) // 2), 
                shield_90_data.dose[Nz // 2 + 25, Ny // 2, Nx//2 + 8:][::2],
                lw=0.0, color='darkgreen', marker='o', 
                markersize=10
            )

            ax.plot(
                linspace(x_min, -1.6, (Nx//2 - 8) // 2), 
                shield_180_data.dose[Nz // 2 + 25, Ny // 2, :Nx//2 - 8][::2],
                lw=0.0, label='180 degrees', 
                color='purple', marker='s', 
                markersize=10
            )
            ax_twin.plot(
                linspace(1.6, x_max, (Nx//2 - 8) // 2), 
                shield_180_data.dose[Nz // 2 + 25, Ny // 2, Nx//2 + 8:][::2],
                lw=0.0, color='purple', marker='s', 
                markersize=10
            )

            ax.plot(
                linspace(x_min, -1.6, (Nx//2 - 8) // 2), 
                shield_270_data.dose[Nz // 2 + 25, Ny // 2, :Nx//2 - 8][::2],
                lw=0.0, label='270 degrees', 
                color='darkorange', marker='^', 
                markersize=10
            )
            ax_twin.plot(
                linspace(1.6, x_max, (Nx//2 - 8) // 2), 
                shield_270_data.dose[Nz // 2 + 25, Ny // 2, Nx//2 + 8:][::2],
                lw=0.0, color='darkorange', marker='^', 
                markersize=10
            )

    ax.legend(loc=0, prop={'size': 18})
    ax.grid(True)
    ax.set_title(title_list[0], fontsize=20)
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

    file_dict = {
        (0, 0): target_dir + '/mlwa_270.phantom_wo_applicator.3ddose',
        (0, 1): target_dir + '/mlwa_270_2.phantom_wo_applicator.3ddose',
        (0, 2): target_dir + '/mlwa_270_3.phantom_wo_applicator.3ddose',
        (0, 3): target_dir + '/mlwa_270_4.phantom_wo_applicator.3ddose',
        (1, 0): target_dir + '/mlwa_270_5.phantom_wo_applicator.3ddose',
        (1, 1): target_dir + '/mlwa_270_6.phantom_wo_applicator.3ddose',
        (1, 2): target_dir + '/mlwa_270_7.phantom_wo_applicator.3ddose',
        (1, 3): target_dir + '/mlwa_270_8.phantom_wo_applicator.3ddose'
    }

    label_dict = {
        (0, 0): 'Base Config.',
        (0, 1): 'Config. 2',
        (0, 2): 'Config. 3',
        (0, 3): 'Config. 4',
        (1, 0): 'Config. 5',
        (1, 1): 'Config. 6',
        (1, 2): 'Config. 7',
        (1, 3): 'Config. 8'
    }

    fig, ax = subplots(
        1, 1, figsize=(10, 10),
        sharex='all', sharey='all'
    )
    fig2, ax2 = subplots(
        1,1, figsize=(10, 10),
        sharex='all', sharey='all'
    )

    for i in xrange(2):
        for j in xrange(4):

            if i == 0 and j == 1:

                full_data = DoseFile(file_dict[(i, j)])

                Nx, Ny, Nz = full_data.shape

                # full_data.dose *= 2.2861e14 # scale to total treatment time
                full_data.dose *= 8.2573429808917e13  # scale to maximum individual dwell time

                full_data.dose /= 5  # normalize to desired dose of 5 Gy
                full_data.dose *= 100  # express in percent. Should see 100% at x=-2cm

                x_min, x_max = full_data.x_extent
                y_min, y_max = full_data.y_extent
                z_min, z_max = full_data.z_extent

                xy_contour = ax.contourf(
                    linspace(x_min, x_max, Nx), linspace(y_min, y_max, Ny),
                    # 25 -> slice 5cm above the origin as that is around where the center of the sources are
                    full_data.dose[Nz // 2 + 25, :, :],
                    # arange(0, 110, 10),
                    cmap=cm.Purples
                )
                # ax.set_title(label_dict[(i ,j)], fontsize=20)

                xz_contour = ax2.contourf(
                    linspace(x_min, x_max, Nx), linspace(z_min, z_max, Nz),
                    full_data.dose[:, Ny // 2, :],
                    # arange(0, 110, 10),
                    cmap=cm.Purples
                )
                # ax2.set_title(label_dict[(i ,j)], fontsize=20)

    ax.grid(True)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.set_xticks(arange(x_min, x_max + 1, 2))
    ax.set_yticks(arange(y_min, y_max + 1, 2))
    
    ax2.xaxis.set_tick_params(labelsize=14)
    ax2.yaxis.set_tick_params(labelsize=14)
    ax2.set_xticks(arange(x_min, x_max + 1, 2))
    ax2.set_yticks(arange(y_min, y_max + 1, 2))
    ax2.grid(True)

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

    fig.savefig(pwd + '/xy_isodose_profile.pdf')

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

    fig2.savefig(pwd + '/xz_isodose_profile.pdf')
    
def isodose_subplots():
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
    target_dir = '/Users/student/research/results_not_to_git' # for work

    file_dict = {
        (0,0):target_dir + '/mlwa_180.phantom_wo_applicator.3ddose',
        (0,1):target_dir + '/mlwa_180_2.phantom_wo_applicator.3ddose',
        (0,2):target_dir + '/mlwa_180_3.phantom_wo_applicator.3ddose',
        (0,3):target_dir + '/mlwa_180_4.phantom_wo_applicator.3ddose',
        (1,0):target_dir + '/mlwa_180_5.phantom_wo_applicator.3ddose',
        (1,1):target_dir + '/mlwa_180_6.phantom_wo_applicator.3ddose',
        (1,2):target_dir + '/mlwa_180_7.phantom_wo_applicator.3ddose',
        (1,3):target_dir + '/mlwa_180_8.phantom_wo_applicator.3ddose'
    }

    label_dict = {
        (0,0):'Base Config.',
        (0,1):'Config. 2',
        (0,2):'Config. 3',
        (0,3):'Config. 4',
        (1,0):'Config. 5',
        (1,1):'Config. 6',
        (1,2):'Config. 7',
        (1,3):'Config. 8'
    }

    fig, ax = subplots(
        2, 4, figsize=(10, 10),
        sharex='all', sharey='all'
    )
    fig2, ax2 = subplots(
        2, 4, figsize=(10, 10),
        sharex='all', sharey='all'
    )

    for i in xrange(2):
        for j in xrange(4):

            full_data = DoseFile(file_dict[(i,j)])

            Nx, Ny, Nz = full_data.shape

            # full_data.dose *= 2.2861e14 # scale to total treatment time
            full_data.dose *= 8.2573429808917e13 # scale to maximum individual dwell time

            full_data.dose /= 5  # normalize to desired dose of 5 Gy
            full_data.dose *= 100  # express in percent. Should see 100% at x=-2cm 

            x_min, x_max = full_data.x_extent
            y_min, y_max = full_data.y_extent
            z_min, z_max = full_data.z_extent

            xy_contour = ax[i,j].contourf(
                linspace(x_min, x_max, Nx), linspace(y_min, y_max, Ny), 
                # 25 -> slice 5cm above the origin as that is around where the center of the sources are in config 2
                full_data.dose[Nz // 2 + 25, :, :],
                arange(0,110,10),
                cmap=cm.Purples
            )

            xz_contour = ax2[i, j].contourf(
                linspace(x_min, x_max, Nx), linspace(z_min, z_max, Nz),
                full_data.dose[:, Ny // 2, :],
                arange(0,110,10),
                cmap=cm.Purples
            )

    for n in xrange(2):
        for m in xrange(4):
            ax[n,m].grid(True)
            ax[n,m].set_title(label_dict[(n,m)],fontsize=14)
            ax[n,m].xaxis.set_tick_params(labelsize=14)
            ax[n,m].yaxis.set_tick_params(labelsize=14)
            ax[n,m].set_xticks(arange(x_min, x_max + 1, 5))
            ax[n,m].set_yticks(arange(y_min, y_max + 1, 5))
            ax2[n,m].grid(True)
            ax2[n,m].set_title(label_dict[(n,m)],fontsize=14)
            ax2[n,m].xaxis.set_tick_params(labelsize=14)
            ax2[n,m].yaxis.set_tick_params(labelsize=14)
            ax2[n,m].set_xticks(arange(x_min, x_max + 1, 5))
            ax2[n,m].set_yticks(arange(y_min, y_max + 1,5))

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

    cax = fig.add_axes([0.91,0.13,0.01,0.7])
    cbar1 = fig.colorbar(
        xy_contour, cax=cax, orientation='vertical',
        ax=ax.ravel().tolist()
    )
    cbar1.set_label('Percentage Isodose (%)',fontsize=24)
    cbar1.ax.tick_params(labelsize=14)
    fig.tight_layout()

    left = 0.1  # the left side of the subplots of the figure
    right = 0.89    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig.savefig(pwd + '/xy_isodose_profile_subplots.pdf')
    
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
    cax2 = fig2.add_axes([0.91,0.13,0.01,0.7])
    cbar2 = fig2.colorbar(
        xz_contour, cax=cax2, orientation='vertical',
        ax=ax.ravel().tolist()
        )
    cbar2.set_label('Percentage Isodose (%)', fontsize=24)
    cbar2.set_clim([0,100])
    cbar2.ax.tick_params(labelsize=14)

    fig2.tight_layout()

    left = 0.1  # the left side of the subplots of the figure
    right = 0.89    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig2.savefig(pwd + '/xz_isodose_profile_subplots.pdf')

if __name__ == "__main__":
    program_name, arguments = argv[0], argv[1:]

    # generate_tdvh()
    # dose_position_plots()
    # rel_dose_position_plot()
    isodose_plot()
    # isodose_subplots()


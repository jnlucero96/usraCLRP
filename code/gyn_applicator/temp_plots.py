#!/usr/bin/env
# author: Joseph Lucero
# created on: 25 May 2018 11:33:46
# purpose: plotting 3ddose files

from __future__ import division

from sys import argv, exit
from os import getcwd

from numpy import (
    linspace, histogram, arange, meshgrid, array, empty, around, sqrt
    )
from scipy.interpolate import RegularGridInterpolator as RGI
from py3ddose import DoseFile, position_to_index

from matplotlib.cm import get_cmap
from matplotlib.style import use
use('seaborn-notebook')
from matplotlib.pyplot import subplot, subplots, figure, close
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
# from mpl_toolkits.mplot3d import Axes3D

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

    shield_type_lst = [
        'Unshielded',
        '90 degree shield',
        '180 degree shield',
        '270 degree shield'
    ]

    vox_size_lst = [
        # '0.5mm',
        '1mm',
        '2mm'
    ]
    vox_size_txt_lst = [
        # '0.5mm',
        '1pt0mm',
        '2pt0mm'
    ]

    title_list = [
        'X - axis',
        'Y - axis',
        'Z - axis'
    ]

    for fig_index in xrange(len(vox_size_lst)):

        fig, ax = subplots(
            3, len(shield_type_lst), figsize=(15, 15),
            sharex='col', sharey='all'
        )

        for index2 in xrange(len(shield_type_lst)):  # iterate through shield types

            full_data = DoseFile(
                file_list[index2].format(vox_size_txt_lst[fig_index])
            )

            Nx, Ny, Nz = full_data.shape

            # scale to maximum individual dwell time
            full_data.dose *= 8.2573429808917e13  

            x_min, x_max = full_data.x_extent
            y_min, y_max = full_data.y_extent
            z_min, z_max = full_data.z_extent

            for index1 in xrange(3):  # iterate through axes

                if index1 == 0:
                    ax[index1, index2].plot(
                        linspace(x_min, x_max,
                                 Nx), full_data.dose[:, Ny // 2, Nz // 2],
                        # yerr=full_data.uncertainty[Nz // 2, Ny // 2, :],
                        lw=3.0
                    )
                elif index1 == 1:
                    ax[index1, index2].plot(
                        linspace(y_min, y_max,
                                 Ny), full_data.dose[Nx // 2, :, Nz // 2],
                        # yerr=full_data.uncertainty[Nz // 2, :, Nx // 2],
                        lw=3.0
                    )
                else:
                    ax[index1, index2].plot(
                        linspace(z_min, z_max,
                                 Nz), full_data.dose[Nx // 2, Ny // 2, :],
                        # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                        lw=3.0
                    )
        for n in xrange(len(title_list)):
            for m in xrange(len(shield_type_lst)):
                ax[n, m].grid(True)
                if n == 0:
                    ax[n, m].set_title(shield_type_lst[m], fontsize=20)
                if m == len(shield_type_lst) - 1:
                    ax[n, m].set_ylabel(title_list[n],fontsize=20)
                    ax[n, m].yaxis.set_label_position("right")
                ax[n, m].xaxis.set_tick_params(labelsize=14)
                ax[n, m].yaxis.set_tick_params(labelsize=14)
                ax[n, m].set_xticks(arange(z_min, z_max + 1, 5))

        fig.text(
            0.01, 0.51, 'Dose (Gy)',
            fontsize=27, rotation='vertical', va='center'
        )
        fig.text(
            0.43, 0.03, 'Position (cm)', fontsize=27, va='center'
        )
        fig.text(
            0.52, 0.95,
            'Absolute Dose vs. Position \n With Volume Correction; ncase = 5E9',
            fontsize=27, va='center', ha='center'
        )
        fig.tight_layout()

        left = 0.125  # the left side of the subplots of the figure
        right = 0.95    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.87     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width reserved for blank space between subplots
        # hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig.savefig(pwd + '/dosage_comparison_' + vox_size_lst[fig_index] + '.pdf')

def dose_inv_position_interpolate_plots(interpolation=False):
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
        '/mlwa_30mmOut_0shield_{0}_sim.phantom_wo_applicator.3ddose',
        target_dir +
        '/mlwa_30mmOut_90shield_{0}_sim.phantom_wo_applicator.3ddose',
        target_dir +
        '/mlwa_30mmOut_180shield_{0}_sim.phantom_wo_applicator.3ddose',
        target_dir +
        '/mlwa_30mmOut_270shield_{0}_sim.phantom_wo_applicator.3ddose'
        # target_dir +
        # '/tg43_{0}_sim.phantom_wo_box.3ddose',
        # target_dir +
        # '/tg43_applicator_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
    ]

    shield_type_lst = [
        'Unshielded',
        '90',
        # '180',
        # '270'
    ]

    vox_size_lst = [
        # '0.5mm',
        '1mm',
        '2mm'
    ]
    vox_size_txt_lst = [
        # '0.5mm',
        '1pt0mm',
        '2pt0mm'
    ]

    title_list = [
        'y-cuts (|x| = {0}, z = {1})',
        'z-cuts (|x| = {0}, y = {1})'
    ]

    for x_pos_desired in [1.62]:
        # for y_pos in around(arange(-2.0, 2.2, 0.2), 1):
        for y_pos_desired in [0.0]:
            for fig_index in xrange(len(vox_size_lst)):
                for index2 in xrange(len(shield_type_lst)):  # iterate through shield types

                    fig = figure(figsize=(15,10))
                    minor_locator = MultipleLocator(0.5)
                    minor_locator2 = MultipleLocator(0.5)
                    gs = GridSpec(ncols=2, nrows=2)

                    ax1 = fig.add_subplot(gs[0, 0])
                    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
                    ax3 = fig.add_subplot(gs[1, :])

                    ax = [ax1, ax2, ax3]

                    full_data = DoseFile(
                        file_list[index2].format(vox_size_txt_lst[fig_index])
                    )

                    # scale to maximum individual dwell time
                    full_data.dose *= 8.2573429808917e13

                    Nx, Ny, Nz = full_data.shape

                    x_min, x_max = full_data.x_extent
                    y_min, y_max = full_data.y_extent
                    z_min, z_max = full_data.z_extent

                    x_pos = array(full_data.positions[0])
                    y_pos = array(full_data.positions[1])
                    z_pos = array(full_data.positions[2])

                    x_pos_mid = (x_pos[1:] + x_pos[:-1]) / 2.0
                    y_pos_mid = (y_pos[1:] + y_pos[:-1]) / 2.0
                    z_pos_mid = (z_pos[1:] + z_pos[:-1]) / 2.0

                    if interpolation:

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

                        z_depths = empty(z_pos.size)

                        for z_depth_index in xrange(z_pos.size):
                            z_depths[z_depth_index] = (
                                dose_matrix[
                                    position_to_index(
                                        x_pos_desired, x_pos),
                                    position_to_index(
                                        y_pos_desired, y_pos),
                                    z_depth_index
                                ] / dose_matrix[
                                    position_to_index(-x_pos_desired,
                                                      x_pos),
                                    position_to_index(
                                        y_pos_desired, y_pos),
                                    z_depth_index
                                ]
                            )

                        ax[0].plot(
                            z_pos,
                            dose_matrix[
                                position_to_index(x_pos_desired, x_pos),
                                # :,
                                position_to_index(y_pos_desired, y_pos),
                                :
                            ],
                            # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                            lw=3.0
                        )
                        ax[1].plot(
                            z_pos,
                            # x_pos, y_pos,
                            dose_matrix[
                                position_to_index(-x_pos_desired, x_pos),
                                # :,
                                position_to_index(y_pos_desired, y_pos),
                                :
                                # z_position_to_index[0.0]
                            ],
                            # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                            lw=3.0
                        )
                        ax[2].plot(
                            z_pos,
                            z_depths * 100,
                            # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                            lw=3.0
                        )

                    else: 
                    
                        dose_matrix = full_data.dose

                        z_depths = empty(z_pos_mid.size)
                        
                        for z_depth_index in xrange(z_pos_mid.size):
                            z_depths[z_depth_index] = (
                                dose_matrix[
                                    position_to_index(x_pos_desired,x_pos_mid),
                                    position_to_index(y_pos_desired,y_pos_mid),
                                    z_depth_index
                                ] / dose_matrix[
                                    position_to_index(-x_pos_desired,x_pos_mid), 
                                    position_to_index(y_pos_desired,y_pos_mid), 
                                    z_depth_index
                                ] 
                            )
                        ax[0].plot(
                            z_pos_mid,
                            dose_matrix[
                                position_to_index(x_pos_desired, x_pos_mid),
                                # :,
                                position_to_index(y_pos_desired, y_pos_mid),
                                :
                            ],
                            # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                            lw=3.0
                        )
                        ax[1].plot(
                            z_pos_mid,
                            # x_pos, y_pos,
                            dose_matrix[
                                position_to_index(-x_pos_desired, x_pos_mid),
                                # :,
                                position_to_index(y_pos_desired, y_pos_mid),
                                :
                                # z_position_to_index[0.0]
                            ],
                            # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                            lw=3.0
                        )
                        ax[2].plot(
                            z_pos_mid,
                            z_depths * 100,
                            # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                            lw=3.0
                        )

                    ax[0].set_ylabel(
                        r"$D_{\mathrm{shielded}}$ (Gy)", fontsize=17
                        )
                    ax[1].set_ylabel(
                        r"$D_{\mathrm{unshielded}}$ (Gy)", fontsize=17
                        )
                    ax[2].set_ylabel(
                        r"$D_{\mathrm{shielded}}\ /\ D_{\mathrm{unshielded}}$ (%)", 
                        fontsize = 17
                        )

                    for m in xrange(3):
                        ax[m].grid(True, which='both')
                        ax[m].set_xticks(arange(z_min, z_max + 1, 1))
                        ax[m].xaxis.set_minor_locator(minor_locator)
                        if m != 2:
                            ax[m].yaxis.set_minor_locator(minor_locator2)
                        ax[m].xaxis.set_tick_params(labelsize=12)
                        ax[m].yaxis.set_tick_params(labelsize=12)
                        ax[m].set_xlim([-10,10])

                    # fig.text(
                    #     0.01, 0.51, 'Ratio at depths',
                    #     fontsize=27, rotation='vertical', va='center'
                    # )
                    fig.text(
                        0.55, 0.03, 'Depth of cut (cm)', fontsize=27, va='center',
                        ha='center'
                    )
                    fig.text(
                        0.52, 0.95,
                        'Dose Ratios vs. Depth \n With Volume Correction; ncase = 1E9',
                        fontsize=27, va='center', ha='center'
                    )
                    fig.tight_layout()

                    left = 0.09  # the left side of the subplots of the figure
                    right = 0.97    # the right side of the subplots of the figure
                    bottom = 0.09   # the bottom of the subplots of the figure
                    top = 0.87     # the top of the subplots of the figure
                    # wspace = 0.2  # the amount of width for blank space between subplots
                    # hspace = 0.2  # the amount of height for white space between subplots

                    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

                    fig.savefig(
                        pwd + '/dosage_inv_comparison_' + vox_size_lst[fig_index] 
                        + '_shield' + str(shield_type_lst[index2])
                        + '_x' + str(x_pos_desired) + '_y' + str(y_pos_desired) 
                        + '.pdf'
                        )

                    close(fig)

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

    title_str = 'Relative Dose vs. Position \n Applicator-as-source; ncase = 5E9'

    fig, ax = subplots(
        1, 1, figsize=(5,3)
        )
    ax_twin = ax.twinx()

    unshielded_full_data = DoseFile(unshielded_file,load_uncertainty=True)
    shield_90_data = DoseFile(shielded_file_90,load_uncertainty=True)
    shield_180_data = DoseFile(shielded_file_180,load_uncertainty=True)
    shield_270_data = DoseFile(shielded_file_270,load_uncertainty=True)

    x_min, x_max = unshielded_full_data.x_extent
    z_min, z_max = unshielded_full_data.z_extent

    x_pos = array(unshielded_full_data.positions[0])
    y_pos = array(unshielded_full_data.positions[1])
    z_pos = array(unshielded_full_data.positions[2])

    x_pos_mid = (x_pos[1:] + x_pos[:-1]) / 2.0
    y_pos_mid = (y_pos[1:] + y_pos[:-1]) / 2.0
    z_pos_mid = (z_pos[1:] + z_pos[:-1]) / 2.0

    # unshielded_interpolated_dose = RGI(
    #     (x_pos_mid, y_pos_mid, z_pos_mid),
    #     unshielded_full_data.dose,
    #     bounds_error=False,
    #     fill_value=None
    # )
    # shield_90_interpolated_dose = RGI(
    #     (x_pos_mid, y_pos_mid, z_pos_mid),
    #     shield_90_data.dose,
    #     bounds_error=False,
    #     fill_value=None
    # )
    # shield_180_interpolated_dose = RGI(
    #     (x_pos_mid, y_pos_mid, z_pos_mid),
    #     shield_180_data.dose,
    #     bounds_error=False,
    #     fill_value=None
    # )
    # shield_270_interpolated_dose = RGI(
    #     (x_pos_mid, y_pos_mid, z_pos_mid),
    #     shield_270_data.dose,
    #     bounds_error=False,
    #     fill_value=None
    # )

    # points = meshgrid(x_pos, y_pos, z_pos)
    # flat = array([m.flatten() for m in points])

    # unshielded_interpolated_dose_matrix = unshielded_interpolated_dose(
    #     flat.transpose()
    # ).reshape(*points[0].shape).transpose((1, 0, 2))
    # shield_90_interpolated_dose_matrix = shield_90_interpolated_dose(
    #     flat.transpose()
    # ).reshape(*points[0].shape).transpose((1, 0, 2))
    # shield_180_interpolated_dose_matrix = shield_180_interpolated_dose(
    #     flat.transpose()
    # ).reshape(*points[0].shape).transpose((1, 0, 2))
    # shield_270_interpolated_dose_matrix = shield_270_interpolated_dose(
    #     flat.transpose()
    # ).reshape(*points[0].shape).transpose((1, 0, 2))

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
        x_pos_mid[::3],
        # shield_90_data.dose[
        shield_90_interpolated_dose_matrix[
            :, 
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        yerr=error_90[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        lw=0.0, label=r'$90^{\circ}$', color='darkgreen',
        markeredgecolor='darkgreen', marker='o',
        markeredgewidth=1, markersize=6, markerfacecolor='None',
        elinewidth=1.0, capsize=0.75
    )
    ax_twin.errorbar(
        x_pos_mid[::3],
        # shield_90_data.dose[
        shield_90_interpolated_dose_matrix[
            :, 
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        yerr=error_90[
            :, 
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        lw=0.0, label=r'$90^{\circ}$', color='darkgreen',
        markeredgecolor='darkgreen', marker='o',
        markeredgewidth=1, markersize=6, markerfacecolor='None',
        elinewidth=1.0, capsize=0.75
    )

    ax.errorbar(
        x_pos_mid[::3],
        # shield_180_data.dose[
        shield_180_interpolated_dose_matrix[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        yerr=error_180[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        lw=0.0, label=r'$180^{\circ}$', color='purple',
        markeredgecolor='purple', marker='s',
        markeredgewidth=1, markersize=6, markerfacecolor='None',
        elinewidth=1.0, capsize=0.75
    )
    ax_twin.errorbar(
        x_pos_mid[::3],
        # shield_180_data.dose[
        shield_180_interpolated_dose_matrix[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        yerr=error_180[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        lw=0.0, label=r'$180^{\circ}$', color='purple',
        markeredgecolor='purple', marker='s',
        markeredgewidth=1, markersize=6, markerfacecolor='None',
        elinewidth=1.0, capsize=0.75
    )

    ax.errorbar(
        x_pos_mid[::3],
        # shield_270_data.dose[
        shield_270_interpolated_dose_matrix[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        yerr=error_270[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        lw=0.0, label=r'$270^{\circ}$', color='darkorange',
        markeredgecolor='darkorange', marker='^',
        markeredgewidth=1, markersize=6, markerfacecolor='None',
        elinewidth=1.0, capsize=0.75
    )
    ax_twin.errorbar(
        x_pos_mid[::3],
        # shield_270_data.dose[
        shield_270_interpolated_dose_matrix[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        yerr=error_270[
            :,
            position_to_index(0.0, y_pos_mid),
            position_to_index(0.0, z_pos_mid)
        ][::3],
        lw=0.0, label=r'$270^{\circ}$', color='darkorange',
        markeredgecolor='darkorange', marker='^',
        markeredgewidth=1, markersize=8, markerfacecolor='None',
        elinewidth=1.0, capsize=0.75
    )

    ax.legend(loc=0, prop={'size': 6})
    # ax.grid(True)
    # ax.set_title('X - axis', fontsize=10)
    ax.set_yticks(arange(0.80, 1.06, 0.05))
    ax.set_ylim([0.80, 1.05])
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.set_xticks(arange(-10, 10 + 1, 1))
    ax.set_xlim([-10, 10])
    ax.axhline(1.0, 0, 0.5, lw=1.5, color='gray')

    ax_twin.set_ylim([0.0, 0.3])
    ax_twin.set_yticks(arange(0.0, 0.31, 0.05))
    ax_twin.yaxis.set_tick_params(labelsize=8)
    ax_twin.vlines([-1.5,1.5],0,1.0)
    ax_twin.fill_between([-1.5, 1.5], 1.05, facecolor='lightgray')

    fig.text(
        0.02, 0.51, 'Relative Dose',
        fontsize=10, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.99, 0.51, 'Effective Transmission',
        fontsize=10, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.5, 0.51, 'A p p l i c a t o r',
        fontsize=11, rotation='vertical', va='center', ha='center'
    )
    fig.text(
        0.52, 0.03, 'Position (cm)', fontsize=10, va='center', ha='center'
    )
    fig.text(
        0.52, 0.93,
        title_str,
        fontsize=10, va='center', ha='center'
    )
    fig.tight_layout()

    left = 0.11  # the left side of the subplots of the figure
    right = 0.89    # the right side of the subplots of the figure
    bottom = 0.15   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    # wspace = 0.2  # the amount of width for blank space between subplots
    # hspace = 0.2  # the amount of height for white space between subplotss

    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig.savefig(pwd + '/relative_dosage_comparison_nb.pdf')
    
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

    title_str = 'Isodose Contours \n Applicator-as-source; ncase = 5E9'

    for fig_index in xrange(len(vox_size_lst)):

        fig, ax = subplots(
            2, 2, figsize=(5, 5),
            sharex='all', sharey='all'
        )
        fig2, ax2 = subplots(
            2, 2, figsize=(5, 5),
            sharex='all', sharey='all'
        )

        for file_index, file in enumerate(file_list):

            full_data = DoseFile(file.format(vox_size_txt_lst[fig_index]))

            y_position_to_index = {
                position:index 
                for index, position in enumerate(full_data.positions[1])
            }
            z_position_to_index = {
                position:index 
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
                full_data.dose[:, :, position_to_index(0.0, z_pos_mid)].transpose(),
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
            ax[ax_x, ax_y].set_title(shield_type_lst[file_index], fontsize=8)

            xz_contour = ax2[ax_x, ax_y].contourf(
                x_pos_mid, z_pos_mid,
                # matplotlib plots column by row (instead of row by column)
                # so transpose data array to account for this
                full_data.dose[:, position_to_index(0.0, y_pos_mid), :].transpose(),
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
            ax2[ax_x, ax_y].set_title(shield_type_lst[file_index], fontsize=8)

        for n in xrange(2):
            for m in xrange(2):

                # ax[n, m].grid(True)
                ax[n, m].xaxis.set_tick_params(
                    labelsize=6, top=True, direction='in'
                    )
                ax[n, m].yaxis.set_tick_params(
                    labelsize=6, right=True, direction='in'
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
                    labelsize=6, top=True, direction='in'
                    )
                ax2[n, m].yaxis.set_tick_params(
                    labelsize=6, right=True, direction='in'
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
            fontsize=10, rotation='vertical', va='center'
        )
        fig.text(
            0.43, 0.03, 'x-axis (cm)', fontsize=10, va='center'
        )
        fig.text(
            0.52, 0.965,
            title_str,
            fontsize=10, va='center', ha='center'
        )

        cax = fig.add_axes([0.88, 0.09, 0.01, 0.79])
        cbar1 = fig.colorbar(
            xy_contour, cax=cax, orientation='vertical',
            ax=ax
        )
        cbar1.set_label('Percentage Isodose (%)', fontsize=10)
        cbar1.ax.tick_params(labelsize=8)
        fig.tight_layout()

        left = 0.1  # the left side of the subplots of the figure
        right = 0.86    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.88     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width reserved for blank space between subplots
        # hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig.savefig(
            pwd + '/xy_isodose_profile_' + vox_size_lst[fig_index] + '_nb.pdf'
        )

        fig2.text(
            0.01, 0.51, 'z-axis (cm)',
            fontsize=10, rotation='vertical', va='center'
        )
        fig2.text(
            0.43, 0.03, 'x-axis (cm)', fontsize=10, va='center'
        )
        fig2.text(
            0.52, 0.965,
            title_str,
            fontsize=10, va='center', ha='center'
        )
        cax2 = fig2.add_axes([0.88, 0.09, 0.01, 0.79])
        cbar2 = fig2.colorbar(
            xz_contour, cax=cax2, orientation='vertical',
            ax=ax
        )
        cbar2.set_label('Percentage Isodose (%)', fontsize=10)
        cbar2.ax.tick_params(labelsize=8)

        fig2.tight_layout()

        left = 0.1  # the left side of the subplots of the figure
        right = 0.86    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.88     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width reserved for blank space between subplots
        # hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig2.savefig(
            pwd + '/xz_isodose_profile_' + vox_size_lst[fig_index] + '_nb.pdf'
        )

def inverse_isodose_plot():
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
        '/mlwa_30mmOut_0shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_30mmOut_90shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_30mmOut_180shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
        target_dir +
        '/mlwa_30mmOut_270shield_{0}_sim.phantom_wo_applicator_wo_box.3ddose',
    ]

    shield_type_lst = [
        'Unshielded',
        r'$90^{\circ}$ shield',
        r'$180^{\circ}$ shield',
        r'$270^{\circ}$ shield'
    ]

    vox_size_lst = [
        '1mm',
        '2mm'
    ]

    vox_size_txt_lst = [
        '1pt0mm',
        '2pt0mm'
    ]

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
                position:index 
                for index, position in enumerate(full_data.positions[1])
            }
            z_position_to_index = {
                position:index 
                for index, position in enumerate(full_data.positions[2])
            }

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

            x_pos = array(full_data.positions[0])[1:] + array(full_data.positions[0])[:-1] / 2.0
            y_pos = array(full_data.positions[1])[1:] + array(full_data.positions[1])[:-1] / 2.0
            z_pos = array(full_data.positions[2])[1:] + array(full_data.positions[2])[:-1] / 2.0

            xy_contour = ax[ax_x, ax_y].contourf(
                # linspace(x_min, x_max, Nx), linspace(y_min, y_max, Ny),
                x_pos, y_pos,
                # matplotlib plots column by row (instead of row by column)
                # so transpose data array to account for this
                full_data.dose[:, :, z_position_to_index[0.0]].transpose(),
                arange(0, 110, 10),
                # [5, 10, 20, 50, 100],
                # cmap=get_cmap('gnuplot')
            )
            ax[ax_x, ax_y].set_title(shield_type_lst[file_index], fontsize=15)

            xz_contour = ax2[ax_x, ax_y].contourf(
                # linspace(x_min, x_max, Nx), linspace(z_min, z_max, Nz),
                x_pos, z_pos,
                # matplotlib plots column by row (instead of row by column)
                # so transpose data array to account for this
                full_data.dose[:, y_position_to_index[0.0], :].transpose(),
                arange(0, 110, 10),
                # [5, 10, 20, 50, 100],
                # cmap=get_cmap('gnuplot')
            )
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
                ax[n, m].set_xticks(arange(x_min, x_max + 1, 2))
                ax[n, m].set_yticks(arange(y_min, y_max + 1, 2))
            
                ax[n, m].vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
                ax[n, m].hlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)

                # ax2[n, m].grid(True)
                ax2[n, m].xaxis.set_tick_params(
                    labelsize=14, top=True, direction='in'
                    )
                ax2[n, m].yaxis.set_tick_params(
                    labelsize=14, right=True, direction='in'
                    )
                ax2[n, m].set_xticks(arange(x_min, x_max + 1, 2))
                ax2[n, m].set_yticks(arange(y_min, y_max + 1, 2))

                ax2[n, m].vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
                

        fig.text(
            0.01, 0.51, 'y-axis (cm)',
            fontsize=27, rotation='vertical', va='center'
        )
        fig.text(
            0.43, 0.03, 'x-axis (cm)', fontsize=27, va='center'
        )
        fig.text(
            0.52, 0.95,
            'Isodose Contours \n With Volume Correction; ncase = 5E9',
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
        top = 0.88     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width reserved for blank space between subplots
        # hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig2.savefig(
            pwd + '/xz_isodose_profile_' + vox_size_lst[fig_index] + '.pdf'
        )

def get_dose_at_points():

    """
    Description:

    Inputs:

    Outputs:
    """

    

if __name__ == "__main__":
    # dose_position_plots()
    # dose_inv_position_interpolate_plots()
    # rel_dose_position_plot()
    isodose_plot()
    # inverse_isodose_plot()

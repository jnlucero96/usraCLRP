#!/usr/bin/env
# author: Joseph Lucero
# created on: 3 May 2018 08:23:46
# purpose: plotting 3ddose files

from __future__ import division

from sys import argv, exit
from os import getcwd

from numpy import linspace, zeros_like, histogram, arange, array
from py3ddose import DoseFile, position_to_index

from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from matplotlib.style import use
use('seaborn-paper')
from matplotlib.pyplot import subplots, close

from normalize import get_conversion_factor

def calc_per_diff(A,B):
    return ((A - B) / B) * 100

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

    shield_type_list = [
        'Unshielded',
        '90 degree shield',
        '180 degree shield',
        '270 degree shield'
    ]

    vxl_size_lst = [
        'Voxel size = (0.5 mm)^{3}',
        'Voxel size = (1.0 mm)^{3}',
        'Voxel size = (2.0 mm)^{3}'
    ]

    fig, ax = subplots(2, 1, figsize=(10, 10), sharex='all', sharey='all')
    fig2, ax2 = subplots(2, 1, figsize=(10, 10), sharex='all', sharey='all')

    for index1 in xrange(2):
        for index2 in xrange(3):

            if index1 == 0:
                pass
            else:

                main_data = DoseFile(file_dict[(index1, index2)])
                # converts to Gy; norm to individual max dwell time
                data_flat = main_data.dose.flatten() * 8.2573429808917e13

                n, bins, __ = ax[index1].hist(
                    data_flat,
                    bins=500,
                    color=color_list[index1],
                    label=shield_type_list[index1],
                    weights=zeros_like(data_flat) + 1. / data_flat.size * 100,
                    alpha=0.4
                )

                n_cum_base = n[::-1].cumsum()[::-1]

                ax2[index1].loglog(
                    bins[:-1], n_cum_base,
                    color=color_list[index1],
                    label=shield_type_list[index1]
                )

    for i in xrange(3):
        ax[i].legend(loc=0, prop={'size': 10})
        ax[i].set_title(vxl_size_lst[i], fontsize=20)
        ax[i].xaxis.set_tick_params(labelsize=17)
        ax[i].yaxis.set_tick_params(labelsize=17)

        ax2[i].legend(loc=0, prop={'size': 10})
        ax2[i].set_title(vxl_size_lst[i], fontsize=20)
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
    # wspace = 0.2  # the amount of width reserved for blank space between subplots
    # hspace = 0.2  # the amount of height reserved for white space between subplotss

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
    # wspace = 0.2  # the amount of width reserved for blank space between subplots
    # hspace = 0.2  # the amount of height reserved for white space between subplotss

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

    file_template = target_dir + \
        '/tg43_applicator_{0}_sim.phantom_wo_applicator_wo_box.3ddose'
    

    vox_size_list = [
        '0.5mm',
        '1mm',
        '2mm'
    ]

    title_list = [
        'X - axis',
        'Y - axis',
        'Z - axis'
    ]

    for fig_index in xrange(3):  # iterate through voxel sizes

        fig, ax = subplots(
            3, 1, figsize=(10, 10),
            sharex='col'
        )

        full_data = DoseFile(
            file_template.format(vox_size_list[fig_index])
        )

        Nx, Ny, Nz = full_data.shape

        # full_data.dose *= 2.2861e14 # scale to total treatment time
        full_data.dose *= 8.2573429808917e13  # scale to maximum individual dwell time

        x_min, x_max = full_data.x_extent
        y_min, y_max = full_data.y_extent
        z_min, z_max = full_data.z_extent

        for index1 in xrange(3):  # iterate through axes

            if index1 == 0:
                ax[index1].plot(
                    linspace(x_min, x_max, Nx), full_data.dose[Nz // 2, Ny // 2, :],
                    # yerr=full_data.uncertainty[Nz // 2, Ny // 2, :],
                    lw=3.0
                )
            elif index1 == 1:
                ax[index1].plot(
                    linspace(y_min, y_max, Ny), full_data.dose[Nz // 2, :, Nx // 2],
                    # yerr=full_data.uncertainty[Nz // 2, :, Nx // 2],
                    lw=3.0
                )
            else:
                ax[index1].plot(
                    linspace(z_min, z_max, Nz), full_data.dose[:, Ny // 2, Nx // 2],
                    # yerr=full_data.uncertainty[:, Ny // 2, Nx // 2],
                    lw=3.0
                )

        for n in xrange(3):
            # ax[n].legend(loc=0,prop={'size':8})
            ax[n].grid(True)
            ax[n].set_title(title_list[n],fontsize=20)
            ax[n].xaxis.set_tick_params(labelsize=14)
            ax[n].yaxis.set_tick_params(labelsize=14)
            # if n == 2:
            #     ax[n].set_xticks(arange(z_min, z_max))

        fig.text(
            0.01, 0.51, 'Dose (Gy)',
            fontsize=27, rotation='vertical', va='center'
        )
        fig.text(
            0.43, 0.03, 'Position (cm)', fontsize=27, va='center'
        )
        fig.text(
            0.52, 0.95,
            'Absolute Dose vs. Position (' + \
            vox_size_list[fig_index] + ')^{3} \n With Volume Correction; ncase = 1E10',
            fontsize=27, va='center', ha='center'
        )
        fig.tight_layout()

        left = 0.125  # the left side of the subplots of the figure
        right = 0.95    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.87     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width reserved for blank space between subplots
        # hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig.subplots_adjust(
            left=left, bottom=bottom, right=right, top=top 
            # wspace=wspace, hspace=hspace
        )

        fig.savefig(
            pwd + '/tg43_applicator_dosage_comparison_' 
            + vox_size_list[fig_index] + '.pdf'
            )


def isodose_plot(place='work', mode='appl'):
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots 
    from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    if place == 'home': 
        target_dir = '/Users/JLucero/MPhysProj/results_not_to_git' # for home
    else:
        target_dir = '/Users/student/research/results_not_to_git'  # for work

    if mode == 'pure':
        file_name_template = '/tg43_{0}pt0mm_sim.phantom_wo_box.3ddose.gz'
    else:
        file_name_template = '/tg43appl_{1}mmOut_{0}pt0mm_sim.phantom_wo_applicator_wo_box.3ddose.gz'

    vox_size_list = [
        '1',
        '2'
    ]

    fig, ax = subplots(
        1, 1, figsize=(10, 10),
        sharex='all', sharey='all'
    )
    fig2, ax2 = subplots(
        1, 1, figsize=(10, 10),
        sharex='all', sharey='all'
    )

    air_kerma = 326.05715 
    air_kerma_per_hist = 1.1584e-13
    max_dwell_time = 0.02917

    for voxel_size in vox_size_list:

        full_data = DoseFile(
            target_dir + file_name_template.format(voxel_size, '30')
        )

        full_data.dose *= get_conversion_factor(
            air_kerma, air_kerma_per_hist, max_dwell_time
            )  # normalize to desired dose of 5 Gy

        full_data.dose /= 5  # normalize to desired dose of 5 Gy
        full_data.dose *= 100  # normalize to desired dose of 5 Gy

        x_min, x_max = full_data.x_extent
        y_min, y_max = full_data.y_extent
        z_min, z_max = full_data.z_extent

        x_pos = array(full_data.positions[0])
        y_pos = array(full_data.positions[1])
        z_pos = array(full_data.positions[2])

        x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2.0
        y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2.0
        z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2.0

        # print x_pos_mid.size, y_pos_mid.size, full_data.dose[:, :, position_to_index(0.0, z_pos)].size
        # exit(0)

        xy_contour = ax.contourf(
            x_pos_mid, y_pos_mid,
            full_data.dose[:, :, position_to_index(0.0, z_pos)].transpose(),
            arange(0, 110, 10),
            # [5, 10, 20, 50, 100]
            # cmap=get_cmap('Purples')
        )
        # ax.set_title(vox_size_list,fontsize=14)

        xz_contour = ax2.contourf(
            x_pos_mid, z_pos_mid,
            full_data.dose[:, position_to_index(0.0, y_pos), :].transpose(),
            arange(0, 110, 10),
            # [5, 10, 20, 50, 100]
            # cmap=get_cmap('Purples')
        )
        # ax2.set_title(vox_size_list, fontsize=14)

        # ax.grid(True)
        ax.xaxis.set_tick_params(labelsize=14)
        ax.yaxis.set_tick_params(labelsize=14)
        ax.set_xticks(arange(-10, 10 + 1, 2))
        ax.set_xlim([-10, 10])
        ax.set_yticks(arange(-10, 10 + 1, 2))
        ax.set_ylim([-10, 10])
        # ax.vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
        # ax.hlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)

        # ax2.grid(True)
        ax2.xaxis.set_tick_params(labelsize=14)
        ax2.yaxis.set_tick_params(labelsize=14)
        ax2.set_xticks(arange(-10, 10 + 1, 2))
        ax2.set_xlim([-10, 10])
        ax2.set_yticks(arange(-10, 10 + 1, 2))
        ax2.set_ylim([-10, 10])
        # ax2.vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)

        fig.text(
            0.01, 0.51, 'y-axis (cm)',
            fontsize=27, rotation='vertical', va='center'
        )
        fig.text(
            0.43, 0.03, 'x-axis (cm)', fontsize=27, va='center'
        )
        # fig.text(
        #     0.52, 0.95,
        #     ,
        #     fontsize=27, va='center', ha='center'
        # )

        cax = fig.add_axes([0.91, 0.09, 0.01, 0.79])
        cbar1 = fig.colorbar(
            xy_contour, cax=cax, orientation='vertical',
            ax=ax, ticks=arange(0, 110, 10)
        )
        cbar1.set_label('Percentage Dose (%)', fontsize=24)

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
            pwd + '/tg43_' + voxel_size + '_xy_isodose_profile_subplots.pdf'
        )

        fig2.text(
            0.01, 0.51, 'z-axis (cm)',
            fontsize=27, rotation='vertical', va='center'
        )
        fig2.text(
            0.43, 0.03, 'x-axis (cm)', fontsize=27, va='center'
        )
        # fig2.text(
        #     0.52, 0.95,
        #     sup_title,
        #     fontsize=27, va='center', ha='center'
        # )
        cax2 = fig2.add_axes([0.90, 0.09, 0.01, 0.86])
        cbar2 = fig2.colorbar(
            xz_contour, cax=cax2, orientation='vertical',
            ax=ax2, ticks=arange(0, 110, 10)
        )
        cbar2.set_label('Percentage Dose (%)', fontsize=24)
        # cbar2.set_clim([0, 100])
        cbar2.ax.tick_params(labelsize=14)

        fig2.tight_layout()

        left = 0.1  # the left side of the subplots of the figure
        right = 0.888    # the right side of the subplots of the figure
        bottom = 0.09   # the bottom of the subplots of the figure
        top = 0.95     # the top of the subplots of the figure
        # wspace = 0.2  # the amount of width reserved for blank space between subplots
        # hspace = 0.2  # the amount of height reserved for white space between subplotss

        fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

        fig2.savefig(
            pwd + '/tg43_' + voxel_size + 'mm_xz_isodose_profile_subplots.pdf'
        )

        close(fig)
        close(fig2)

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

    # file_name_template = target_dir + \
    #     '/tg43_applicator_{0}_sim.phantom_wo_applicator_wo_box.3ddose'
    file_name_template = target_dir + \
        '/tg43_{0}_{1}.phantom_wo_box.3ddose'

    vox_size_list = [
        # '0.5mm',
        # '1mm',
        '2mm'
    ]
    vox_size_txt_list = [
        # '0.5mm',
        # '1pt0mm',
        # '2pt0mm'
        '2mm'
    ]

    sources = [
        'Flexisource',
        'microSelectronV1',
        'microSelectronV2'
    ]

    fig, ax = subplots(
        1, 1, figsize=(10, 10),
        sharex='all', sharey='all'
    )
    fig2, ax2 = subplots(
        1, 1, figsize=(10, 10),
        sharex='all', sharey='all'
    )
    fig3, ax3 = subplots(
        1, 1, figsize=(10, 10),
        sharex='all', sharey='all'
    )

    # source = 'Flexi-V2Compare'
    # source = 'Flexi-V1Compare'
    source = 'V2-V1Compare'
    
    if source is 'V2-V1Compare':
        print "Source is:", source
        data1 = DoseFile(
            target_dir + '/tg43_microSelectronV2_2mm.phantom_wo_box.3ddose.gz'
        )
        data2 = DoseFile(
            target_dir + '/tg43_microSelectronV1_2mm.phantom_wo_box.3ddose.gz'
        )
        comparison_matrix = calc_per_diff(
            data1.dose * 82583025662064.77,
            data2.dose * 82105378673169.89
            )
        sup_title = 'microSelectron-v1 to microSelectron-v2\n Dose Comparison'
    elif source is 'Flexi-V1Compare':
        print "Source is:", source
        data1 = DoseFile(
            target_dir + '/tg43_Flexisource_2mm.phantom_wo_box.3ddose.gz'
        )
        data2 = DoseFile(
            target_dir + '/tg43_microSelectronV1_2mm.phantom_wo_box.3ddose.gz'
        )
        comparison_matrix = calc_per_diff(
            data1.dose * 79906971237618.34,
            data2.dose * 82105378673169.89
            )
        sup_title = 'Flexisource to microSelectron-v1\n Dose Comparison'
    elif source is 'Flexi-V2Compare':
        print "Source is:", source
        data1 = DoseFile(
            target_dir + '/tg43_Flexisource_2mm.phantom_wo_box.3ddose.gz'
        )
        data2 = DoseFile(
            target_dir + '/tg43_microSelectronV2_2mm.phantom_wo_box.3ddose.gz'
        )
        comparison_matrix = calc_per_diff(
            data1.dose * 79906971237618.34,
            data2.dose * 82583025662064.77
            )
        sup_title = 'Flexisource to microSelectron-v2\n Dose Comparison'
    else:
        print "Prompt not understood."; exit(1)
    
    Nx, Ny, Nz = data1.shape

    x_min, x_max = data1.x_extent
    y_min, y_max = data1.y_extent
    z_min, z_max = data1.z_extent

    x_pos = array(data1.positions[0])
    y_pos = array(data1.positions[1])
    z_pos = array(data1.positions[2])

    x_pos_mid = (x_pos[:-1] + x_pos[1:]) / 2.0
    y_pos_mid = (y_pos[:-1] + y_pos[1:]) / 2.0
    z_pos_mid = (z_pos[:-1] + z_pos[1:]) / 2.0

    xy_contour = ax.contourf(
        x_pos_mid, y_pos_mid, 
        comparison_matrix[:,:,position_to_index(0.0,z_pos)].transpose(),
        arange(-7, 10.5, 0.5),
        # [5, 10, 20, 50, 100]
        cmap=get_cmap('spectral')
    )
    # ax.set_title(vox_size_list,fontsize=14)

    xz_contour = ax2.contourf(
        x_pos_mid, z_pos_mid, 
        comparison_matrix[:, position_to_index(0.0, y_pos), :].transpose(),
        arange(-7, 10.5, 0.5),
        # [5, 10, 20, 50, 100]
        cmap=get_cmap('spectral')
    )
    # ax2.set_title(vox_size_list, fontsize=14)

    xy_hist_data = ax3.hist(
        comparison_matrix[:, :, :].flatten(),
        bins='auto', density=True
    )

    # ax.grid(True)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.set_xticks(arange(-10, 10 + 1, 2))
    ax.set_xlim([-10, 10])
    ax.set_yticks(arange(-10, 10 + 1, 2))
    ax.set_ylim([-10, 10])
    # ax.vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)
    # ax.hlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)

    # ax2.grid(True)
    ax2.xaxis.set_tick_params(labelsize=14)
    ax2.yaxis.set_tick_params(labelsize=14)
    ax2.set_xticks(arange(-10, 10 + 1, 2))
    ax2.set_xlim([-10, 10])
    ax2.set_yticks(arange(-10, 10 + 1, 2))
    ax2.set_ylim([-10, 10])
    # ax2.vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)

    ax3.grid(True)
    ax3.xaxis.set_tick_params(labelsize=14)
    ax3.yaxis.set_tick_params(labelsize=14)
    ax3.set_xticks(arange(-5, 5 + 1, 1))
    ax3.set_xlim([-5, 5])
    ax3.set_yticks(arange(0, 2.1 + 1, 0.5))
    ax3.set_ylim([0, 1])
    # ax2.vlines([-2, 2], -10, 10, linestyles='dashed', lw=2.0)

    fig.text(
        0.01, 0.51, 'y-axis (cm)',
        fontsize=27, rotation='vertical', va='center'
    )
    fig.text(
        0.43, 0.03, 'x-axis (cm)', fontsize=27, va='center'
    )
    fig.text(
        0.52, 0.95,
        sup_title,
        fontsize=27, va='center', ha='center'
    )

    cax = fig.add_axes([0.91, 0.09, 0.01, 0.79])
    cbar1 = fig.colorbar(
        xy_contour, cax=cax, orientation='vertical',
        ax=ax, ticks=arange(0,15,1)
    )
    cbar1.set_label('Percentage difference (%)', fontsize=24)
    
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
        pwd + '/tg43_' + source +'_xy_isodose_profile_subplots.pdf'
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
        sup_title,
        fontsize=27, va='center', ha='center'
    )
    cax2 = fig2.add_axes([0.90, 0.09, 0.01, 0.79])
    cbar2 = fig2.colorbar(
        xz_contour, cax=cax2, orientation='vertical',
        ax=ax2, ticks=arange(0,15,1)
    )
    cbar2.set_label('Percentage difference (%)', fontsize=24)
    # cbar2.set_clim([0, 100])
    cbar2.ax.tick_params(labelsize=14)

    fig2.tight_layout()

    left = 0.1  # the left side of the subplots of the figure
    right = 0.888    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.88     # the top of the subplots of the figure
    # wspace = 0.2  # the amount of width reserved for blank space between subplots
    # hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig2.savefig(
        pwd + '/tg43_' + source + '_xz_isodose_profile_subplots.pdf'
        )

    fig3.text(
        0.03, 0.51, 'Probability density',
        fontsize=27, rotation='vertical', va='center', ha='center'
    )
    fig3.text(
        0.56, 0.03, 'Percentage difference (%)', fontsize=27, va='center',
        ha='center'
    )
    fig3.text(
        0.52, 0.95,
        sup_title,
        fontsize=27, va='center', ha='center'
    )

    fig3.tight_layout()

    left = 0.1  # the left side of the subplots of the figure
    right = 0.95    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.88     # the top of the subplots of the figure
    # wspace = 0.2  # the amount of width reserved for blank space between subplots
    # hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig3.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig3.savefig(
        pwd + '/tg43_' + source + '_comparison_histograms.pdf'
    )
    
if __name__ == "__main__":
    program_name, arguments = argv[0], argv[1:]

    # generate_tdvh_mlwa()
    # generate_tdvh_tg43()
    # dose_position_plots()
    # isodose_plot()
    isodose_plot_compare()


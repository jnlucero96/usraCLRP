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

def generate_plot(file_list):

    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list

    """

    pwd = getcwd()

    file_name_1, file_name_2 = file_list[0], file_list[1]

    full_data = DoseFile(file_name_1,load_uncertainty=True)
    phantom_seed_data = DoseFile(file_name_2,load_uncertainty=True)

    fig_x, ax_x = subplots(
        20, 20, figsize=(30, 30), 
        sharex='all', sharey='all'
        )

    x_min, x_max = full_data.x_extent

    for j in xrange(20):
        for k in xrange(20):

            Nx = full_data.dose[:, j, k].size

            ax_x[j, k].plot(linspace(x_min, x_max, Nx),
                            full_data.dose[:, j, k], lw=3.0, label='Full Geometry')
            ax_x[j, k].plot(linspace(x_min, x_max, Nx),
                            phantom_seed_data.dose[:, j, k], lw=3.0, ls='--', label='Scored Phantom only')
            # ax_y.legend(loc=0,prop={'size':15})
            # ax_y.set_ylabel('Dose (Gy)',fontsize=20)
            # ax_y.set_xlabel('Y Axis distance (cm)',fontsize=20)
            # ax_y.set_title('Comparison (Full vs Scoring Geometry)',fontsize=20)
            ax_x[j, k].set_xticks([])
            ax_x[j, k].set_yticks([])

    fig_x.tight_layout()
    fig_x.savefig(pwd + '/x_dosage_comparison.pdf')

    fig_y, ax_y = subplots(20,20,figsize=(30,30),sharex='all',sharey='all')

    y_min, y_max = full_data.y_extent

    for i in xrange(20):
        for k in xrange(20):

            Ny = full_data.dose[i,:,k].size

            ax_y[i,k].plot(linspace(y_min, y_max, Ny), full_data.dose[i,:,k], lw=3.0, label='Full Geometry')
            ax_y[i,k].plot(linspace(y_min, y_max, Ny),
                    phantom_seed_data.dose[i,:,k], lw=3.0, ls='--', label='Scored Phantom only')
            # ax_y.legend(loc=0,prop={'size':15})
            # ax_y.set_ylabel('Dose (Gy)',fontsize=20)
            # ax_y.set_xlabel('Y Axis distance (cm)',fontsize=20)
            # ax_y.set_title('Comparison (Full vs Scoring Geometry)',fontsize=20)
            ax_y[j, k].set_xticks([])
            ax_y[j, k].set_yticks([])

    fig_y.tight_layout()
    fig_y.savefig(pwd + '/y_dosage_comparison.pdf')

    fig_z, ax_z = subplots(20, 20, figsize=(
        30, 30), sharex='all', sharey='all')

    z_min, z_max = full_data.z_extent

    for i in xrange(20):
        for j in xrange(20):

            Nz = full_data.dose[i, j, :].size

            ax_z[i, j].plot(linspace(z_min, z_max, Nz),
                            full_data.dose[i, j, :], lw=3.0, label='Full Geometry')
            ax_z[i, j].plot(linspace(z_min, z_max, Nz),
                            phantom_seed_data.dose[i, j, :], lw=3.0, ls='--', label='Scored Phantom only')
            # ax_y.legend(loc=0,prop={'size':15})
            # ax_y.set_ylabel('Dose (Gy)',fontsize=20)
            # ax_y.set_xlabel('Y Axis distance (cm)',fontsize=20)
            # ax_y.set_title('Comparison (Full vs Scoring Geometry)',fontsize=20)
            ax_z[i, j].set_xticks([])
            ax_z[i, j].set_yticks([])

    fig_z.tight_layout()
    fig_z.savefig(pwd + '/z_dosage_comparison.pdf')

def generate_plot2(file_list):

    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list

    """

    pwd = getcwd()

    file_name_1, file_name_2, file_name_3 = file_list[0], file_list[1], file_list[2]

    full_prostate_data = DoseFile(file_name_1,load_uncertainty=True)
    phantom_seed_data = DoseFile(file_name_2,load_uncertainty=True)
    water_pad_data = DoseFile(file_name_3, load_uncertainty=True)

    # full_set = hstack(
    #     (
    #         full_prostate_data.dose[:, :, :].flatten(), 
    #         phantom_seed_data.dose[:, :, :].flatten(), 
    #         water_pad_data.dose[:, :, :].flatten()
    #         )
    #     )

    # hist(full_set,bins='auto')
    
    fig_x, ax_x = subplots(
        full_prostate_data.dose[0, :, 0].size, full_prostate_data.dose[0, 0, :].size, 
        figsize=(30, 30), sharex='all', sharey='all'
        )

    x_min, x_max = full_prostate_data.x_extent

    for j in xrange(full_prostate_data.dose[0, :, 0].size):
        for k in xrange(full_prostate_data.dose[0, 0, :].size):

            Nx = full_prostate_data.dose[:, j, k].size

            ax_x[j, k].plot(linspace(x_min, x_max, Nx),
                            full_prostate_data.dose[:, j, k], lw=3.0, label='Full Prostate Geometry')
            ax_x[j, k].plot(linspace(x_min, x_max, Nx),
                            phantom_seed_data.dose[:, j, k], lw=3.0, ls='--', label='Scored Phantom only')
            ax_x[j, k].plot(linspace(x_min, x_max, Nx),
                            water_pad_data.dose[:, j, k], lw=3.0, ls='dotted', label='Water Padded Data')
            # ax_y.legend(loc=0,prop={'size':15})
            # ax_y.set_ylabel('Dose (Gy)',fontsize=20)
            # ax_y.set_xlabel('Y Axis distance (cm)',fontsize=20)
            # ax_y.set_title('Comparison (Full vs Scoring Geometry)',fontsize=20)
            ax_x[j, k].set_xticks([])
            ax_x[j, k].set_yticks([])

    fig_x.tight_layout()
    fig_x.savefig(pwd + '/x_dosage_comparison.pdf')

    fig_y, ax_y = subplots(
        full_prostate_data.dose[:, 0, 0].size, full_prostate_data.dose[0, 0, :].size,
        figsize=(30, 30), sharex='all', sharey='all'
        )

    y_min, y_max = full_prostate_data.y_extent

    for i in xrange(full_prostate_data.dose[:, 0, 0].size):
        for k in xrange(full_prostate_data.dose[0, 0, :].size):

            Ny = full_prostate_data.dose[i,:,k].size

            ax_y[i,k].plot(linspace(y_min, y_max, Ny), full_prostate_data.dose[i,:,k], lw=3.0, label='Full Prostate Geometry')
            ax_y[i,k].plot(linspace(y_min, y_max, Ny),
                    phantom_seed_data.dose[i,:,k], lw=3.0, ls='--', label='Scored Phantom only')
            ax_y[i,k].plot(linspace(y_min, y_max, Ny),
                    water_pad_data.dose[i,:,k], lw=3.0, ls='dotted', label='Water Padded Data')
            # ax_y.legend(loc=0,prop={'size':15})
            # ax_y.set_ylabel('Dose (Gy)',fontsize=20)
            # ax_y.set_xlabel('Y Axis distance (cm)',fontsize=20)
            # ax_y.set_title('Comparison (Full vs Scoring Geometry)',fontsize=20)
            ax_y[j, k].set_xticks([])
            ax_y[j, k].set_yticks([])

    fig_y.tight_layout()
    fig_y.savefig(pwd + '/y_dosage_comparison.pdf')

    fig_z, ax_z = subplots(
        full_prostate_data.dose[:, 0, 0].size, full_prostate_data.dose[0, :, 0].size,
        figsize=(30, 30), sharex='all', sharey='all'
        )

    z_min, z_max = full_prostate_data.z_extent

    for i in xrange(full_prostate_data.dose[:, 0, 0].size):
        for j in xrange(full_prostate_data.dose[0, :, 0].size):

            Nz = full_prostate_data.dose[i, j, :].size

            ax_z[i, j].plot(linspace(z_min, z_max, Nz),
                            full_prostate_data.dose[i, j, :], lw=3.0, label='Full Prostate Geometry')
            ax_z[i, j].plot(linspace(z_min, z_max, Nz),
                            phantom_seed_data.dose[i, j, :], lw=3.0, ls='--', label='Scored Phantom only')
            ax_z[i, j].plot(linspace(z_min, z_max, Nz),
                            water_pad_data.dose[i, j, :], lw=3.0, ls='dotted', label='Water padded data')
            # ax_y.legend(loc=0,prop={'size':15})
            # ax_y.set_ylabel('Dose (Gy)',fontsize=20)
            # ax_y.set_xlabel('Y Axis distance (cm)',fontsize=20)
            # ax_y.set_title('Comparison (Full vs Scoring Geometry)',fontsize=20)
            ax_z[i,j].set_xticks([])
            ax_z[i,j].set_yticks([])

    fig_z.tight_layout()
    fig_z.savefig(pwd + '/z_dosage_comparison.pdf')

def generate_dose_histogram(file_list):

    pwd = getcwd()

    file_name_1, file_name_2, file_name_3 = file_list[0], file_list[1], file_list[2]

    full_prostate_data = DoseFile(file_name_1, load_uncertainty=True)
    phantom_seed_data = DoseFile(file_name_2, load_uncertainty=True)
    water_pad_data = DoseFile(file_name_3, load_uncertainty=True)

    fig,ax = subplots(4,1,figsize=(10,10),sharex='col',sharey='all')

    __, bins_full, __ = ax[0].hist(
        full_prostate_data.dose[:, :, :].flatten(),
        bins='auto', color='darkgreen',label='All Prostate'
    )
    ax[0].legend(loc=0,prop={'size':12})


    __, bins_phantom, __ = ax[1].hist(
        phantom_seed_data.dose[:, :, :].flatten(),
        bins='auto', color='purple',label='Scoring Phantom Only'
    )
    ax[1].legend(loc=0, prop={'size': 12})

    __, bins_water, __ = ax[2].hist(
        water_pad_data.dose[:, :, :].flatten(),
        bins='auto', color='darkorange',label='Water Padding'
    )
    ax[2].legend(loc=0, prop={'size': 12})
    
    ax[3].hist(
        full_prostate_data.dose[:, :, :].flatten(),
        color='darkgreen',
        bins='auto',
        alpha=0.3,
        label='All Prostate'
    )
    ax[3].hist(
        phantom_seed_data.dose[:, :, :].flatten(),
        color='purple',
        bins='auto',
        alpha=0.3,
        label='Scoring Phantom Only'
    )
    ax[3].hist(
        water_pad_data.dose[:, :, :].flatten(),
        color='darkorange',
        bins='auto',
        alpha=0.3,
        label='Water Padded'
    )
    ax[3].legend(loc=0,prop={'size':13})

    for i in xrange(4):
        ax[i].grid(True)

    fig.text(
        0.03,0.51,'Volume (number of voxels)',fontsize=20,rotation='vertical',va='center'
    )
    fig.text(
        0.45,0.03,'Dose (Gy)',fontsize=20,va='center'
    )
    fig.text(
        0.425,0.97,'Dose Histograms',fontsize=20,va='center'
    )

    fig.tight_layout()

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.935      # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig.subplots_adjust(left=left, top=top, right=right, bottom=bottom)
    fig.savefig(pwd + '/dose_histogram.pdf')

    fig2,ax2 = subplots(4,1,figsize=(10,10),sharex='col',sharey='all')

    ax2[0].hist(
        full_prostate_data.dose[:, :, :].flatten(),
        bins='auto', color='darkgreen',cumulative=-1, label='All Prostate'
    )
    ax2[0].legend(loc=0, prop={'size': 12})

    ax2[1].hist(
        phantom_seed_data.dose[:, :, :].flatten(),
        bins='auto', color='purple', cumulative=-1, label='Scoring Phantom Only'
    )
    ax2[1].legend(loc=0, prop={'size': 12})

    ax2[2].hist(
        water_pad_data.dose[:, :, :].flatten(),
        bins='auto', color='darkorange', cumulative=-1, label='Water Padded'
    )
    ax2[2].legend(loc=0, prop={'size': 12})

    ax2[3].hist(
        full_prostate_data.dose[:, :, :].flatten(),
        cumulative=-1,
        color='darkgreen',
        bins='auto',
        alpha=0.3,
        label='All Prostate'
    )
    ax2[3].hist(
        phantom_seed_data.dose[:, :, :].flatten(),
        cumulative=-1,
        color='purple',
        bins='auto',
        alpha=0.3,
        label='Scoring Phantom Only'
    )
    ax2[3].hist(
        water_pad_data.dose[:, :, :].flatten(),
        cumulative=-1,
        color='darkorange',
        bins='auto',
        alpha=0.3,
        label='Water Padded'
    )
    ax2[3].legend(loc=0,prop={'size':13})

    fig2.text(
        0.03, 0.51, 'Volume (number of voxels)', 
        fontsize=20, rotation='vertical', va='center'
    )
    fig2.text(
        0.45, 0.03, 'Dose (Gy)', 
        fontsize=20, va='center'
    )
    fig2.text(
        0.35, 0.97, 'Cumulative Dose Histograms', 
        fontsize=20, va='center'
    )

    for i in xrange(4):
        ax2[i].grid(True)

    fig.tight_layout()

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.935      # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig2.subplots_adjust(left=left, top=top, right=right, bottom=bottom)
    fig2.savefig(pwd + '/cum_dose_histogram.pdf')

def generate_tdvh():

    pwd = getcwd()

    target_dir = '/Users/JLucero/MPhysProj/results_not_to_git'

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
        data_base_flat = main_base_data.dose.flatten() * 2.2681e14 # converts to Gy
        
        main_90_data = DoseFile(file_90_dict[index])
        data_90_flat = main_90_data.dose.flatten() * 2.2681e14 # converts to Gy
        
        main_180_data = DoseFile(file_180_dict[index])
        data_180_flat = main_180_data.dose.flatten() * 2.2681e14 # converts to Gy

        main_270_data = DoseFile(file_270_dict[index])
        data_270_flat = main_270_data.dose.flatten() * 2.2681e14  # converts to Gy

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

    target_dir = '/Users/JLucero/MPhysProj/results_not_to_git'

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

        full_data.dose *= 2.2681e14

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
    
def isodose_plots():
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots 
    from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list
    """

    pwd = getcwd()

    target_dir = '/Users/JLucero/MPhysProj/results_not_to_git'

    file_dict = {
        (0,0):target_dir + '/mlwa_270.phantom_wo_applicator.3ddose',
        (0,1):target_dir + '/mlwa_270_2.phantom_wo_applicator.3ddose',
        (0,2):target_dir + '/mlwa_270_3.phantom_wo_applicator.3ddose',
        (0,3):target_dir + '/mlwa_270_4.phantom_wo_applicator.3ddose',
        (1,0):target_dir + '/mlwa_270_5.phantom_wo_applicator.3ddose',
        (1,1):target_dir + '/mlwa_270_6.phantom_wo_applicator.3ddose',
        (1,2):target_dir + '/mlwa_270_7.phantom_wo_applicator.3ddose',
        (1,3):target_dir + '/mlwa_270_8.phantom_wo_applicator.3ddose'
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
        2, 4, figsize=(10, 10),
        sharex='all',sharey='all'
    )
    fig2, ax2 = subplots(
        2, 4, figsize=(10, 10),
        sharex='col'
    )

    for i in xrange(2):
        for j in xrange(4):

            full_data = DoseFile(file_dict[(i,j)])

            Nx, Ny, Nz = full_data.shape

            full_data.dose *= 2.2681e14

            x_min, x_max = full_data.x_extent
            y_min, y_max = full_data.y_extent
            z_min, z_max = full_data.z_extent

            # plot_extent = [x_min, x_max, y_min, y_max]

            xy_contour = ax[i,j].contourf(
                linspace(x_min, x_max, Nx), linspace(y_min, y_max, Ny), 
                full_data.dose[Nz // 2, :, :],
                arange(0,4.1,0.1),
                cmap=cm.Purples
                # colors='purple'
            )
            # ax[i,j].clabel(xy_contour,inline=1)

            xz_contour = ax2[i, j].contourf(
                linspace(x_min, x_max, Nx), linspace(z_min, z_max, Nz),
                full_data.dose[:, Ny // 2, :],
                arange(0,4.1,0.1),
                cmap=cm.Purples
                # colors='purple'
            )
            # ax2[i,j].clabel(xz_contour,inline=1)
            


    for n in xrange(2):
        for m in xrange(4):
            ax[n,m].grid(True)
            ax[n,m].set_title(label_dict[(n,m)],fontsize=14)
            ax[n,m].xaxis.set_tick_params(labelsize=10)
            ax[n,m].yaxis.set_tick_params(labelsize=10)
            # ax[n,m].set_xticks(arange(x_min,x_max))
            # ax[n,m].set_yticks(arange(y_min,y_max))
            ax2[n,m].grid(True)
            ax2[n,m].set_title(label_dict[(n,m)],fontsize=14)
            ax2[n,m].xaxis.set_tick_params(labelsize=10)
            ax2[n,m].yaxis.set_tick_params(labelsize=10)
            # ax[n,m].set_xticks(arange(x_min,x_max))
            # ax[n,m].set_yticks(arange(y_min,y_max))

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

    cax = fig.add_axes([0.95,0.13,0.01,0.7])
    fig.colorbar(
        xy_contour, cax=cax, orientation='vertical',
        ax=ax.ravel().tolist()
        )
    fig.tight_layout()

    left = 0.1  # the left side of the subplots of the figure
    right = 0.91    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig.savefig(pwd + '/xy_isodose_profiles.pdf')
    
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
    cax2 = fig2.add_axes([0.95,0.13,0.01,0.7])
    fig2.colorbar(xz_contour,cax=cax2,orientation='vertical')

    fig2.tight_layout()

    left = 0.1  # the left side of the subplots of the figure
    right = 0.91    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.87     # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig2.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    fig2.savefig(pwd + '/xz_isodose_profiles.pdf')

if __name__ == "__main__":
    program_name, arguments = argv[0], argv[1:]

    # generate_plot(arguments)
    # generate_plot2(arguments)
    # dose_position_plots()
    # generate_dose_histogram(arguments)
    # generate_tdvh()
    isodose_plots()


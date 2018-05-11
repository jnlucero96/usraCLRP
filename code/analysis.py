#!/usr/bin/env
# author: Joseph Lucero
# created on: 3 May 2018 08:23:46
# purpose: plotting 3ddose files

from __future__ import division

from sys import argv, exit
from os import getcwd

from numpy import linspace, zeros_like, histogram
from py3ddose import DoseFile

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

    data_file_name = (
        pwd + '/m_lymper_source_with_appl_2.phantom_wo_applicator.3ddose'
    )

    main_data = DoseFile(data_file_name, load_uncertainty=True)

    data_flat = main_data.dose[:, :, :].flatten() * 2.6181e14 # converts to Gy

    fig, ax = subplots(2, 1, figsize=(10, 10))

    x_num, y_num, z_num = main_data.shape

    n, bins, __ = ax[0].hist(
        data_flat,
        bins=500,
        color='darkgreen',
        label='Absolute DVH',
        weights=zeros_like(data_flat) + 1. / data_flat.size * 100
        )
    ax[0].legend(loc=0, prop={'size': 12})
    ax[0].grid(True)
    
    n_cum, bins_cum, __ = ax[1].hist(
        data_flat,
        bins=500,
        color='darkgreen',
        label='Cumulative DVH',
        weights=zeros_like(data_flat) + 1. / data_flat.size * 100,
        cumulative=-1
        )
    ax[1].legend(loc=0, prop={'size': 12})
    ax[1].grid(True)

    fig.text(
        0.03, 0.51, 'Volume (%)', 
        fontsize=22, rotation='vertical', va='center'
    )
    fig.text(
        0.45, 0.03, 'Dose (Gy)', fontsize=22, va='center'
    )
    fig.text(
        0.52, 0.95, 
        'Dose Volume Histograms \n With Volume Correction; ncase = 1E8',
        fontsize=22, va='center', ha='center'
    )

    fig.tight_layout()

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.91     # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig.subplots_adjust(left=left, top=top, right=right, bottom=bottom)
    fig.savefig(pwd + '/dose_volume_histogram.pdf')

    fig2, ax2 = subplots(2, 1, figsize=(10,10))

    ax2[0].loglog(
        bins[:-1], n, 
        color='darkgreen',
        label='Absolute DVH'
        )
    ax2[0].legend(loc=0,prop={'size':12})
    
    ax2[1].loglog(
        bins_cum[:-1], n_cum, 
        color='darkgreen',
        label='Absolute DVH'
        )
    ax2[1].legend(loc=0,prop={'size':12})

    fig2.text(
        0.03, 0.51, 'Volume (%)',
        fontsize=22, rotation='vertical', va='center'
    )
    fig2.text(
        0.45, 0.03, 'Dose (Gy)', fontsize=22, va='center'
    )
    fig2.text(
        0.52, 0.95,
        'Dose Volume Histograms \n With Volume Correction; ncase = 1E8',
        fontsize=22, va='center', ha='center'
    )

    fig2.tight_layout()

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.09   # the bottom of the subplots of the figure
    top = 0.91     # the top of the subplots of the figure
    wspace = 0.2  # the amount of width reserved for blank space between subplots
    hspace = 0.2  # the amount of height reserved for white space between subplotss

    fig2.subplots_adjust(left=left, top=top, right=right, bottom=bottom)
    fig2.savefig(pwd + '/dose_volume_histogram_loglog.pdf')


def generate_plot3(file_list):
    """
    Description:
    Takes any number of .3ddose files and plots a plethora of diagnostic plots from the data

    Inputs:
    :name list_file: a list of file names that are to be loaded
    :type list_file: list

    """

    pwd = getcwd()

    file_name_1 = file_list[0]

    full_data = DoseFile(file_name_1, load_uncertainty=True)

    Nx, Ny, Nz = full_data.shape

    fig_x, ax_x = subplots(
        Ny, Nz, figsize=(30, 30),
        sharex='all', sharey='all'
    )

    x_min, x_max = full_data.x_extent

    # for j in xrange(Ny):
    #     for k in xrange(Nz):

    ax_x[j, k].plot(linspace(x_min, x_max, Nx),
                    full_data.dose[:, 0, 0], lw=3.0)
    # ax_x[j, k].set_xticks([])
    # ax_x[j, k].set_yticks([])

    fig_x.tight_layout()
    fig_x.savefig(pwd + '/x_dosage_comparison.pdf')

    fig_y, ax_y = subplots(
        Nx, Nz, 
        figsize=(30, 30), 
        sharex='all', sharey='all'
        )

    y_min, y_max = full_data.y_extent

    # for i in xrange(Nx):
    #     for k in xrange(Nz):

    ax_y[i, k].plot(linspace(y_min, y_max, Ny),
                    full_data.dose[0, :, 0], lw=3.0)
    # ax_y[j, k].set_xticks([])
    # ax_y[j, k].set_yticks([])

    fig_y.tight_layout()
    fig_y.savefig(pwd + '/y_dosage_comparison.pdf')

    fig_z, ax_z = subplots(
        Nx, Ny, 
        figsize=(30, 30),
        sharex='all', sharey='all'
        )

    z_min, z_max = full_data.z_extent

    # for i in xrange(Nx):
    #     for j in xrange(Ny):

    ax_z[i, j].plot(linspace(z_min, z_max, Nz),
                    full_data.dose[0, 0, :], lw=3.0)
    # ax_z[i, j].set_xticks([])
    # ax_z[i, j].set_yticks([])

    fig_z.tight_layout()
    fig_z.savefig(pwd + '/z_dosage_comparison.pdf')

if __name__ == "__main__":
    program_name, arguments = argv[0], argv[1:]

    # generate_plot(arguments)
    # generate_plot2(arguments)
    generate_plot3(arguments)
    # generate_dose_histogram(arguments)
    # generate_tdvh()


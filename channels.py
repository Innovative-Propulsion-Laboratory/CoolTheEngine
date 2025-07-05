# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 21:46:19 2020

@author: julien
"""
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import cte_tools as t
from scipy.interpolate import interp1d


def canaux(profile_data, width_data, height_data, angle_data, wall_thickness,
           nbc, x_chamber_throat_exit, plot_detail, write_in_csv, figure_dpi, plot_dir=None):
    """
    This function computes the caracteristics of channels on each point
    by interpolation between given values at injection plate (inj), end of cylindrical chamber (conv), 
    throat (col) and extremity of the nozzle (div).
    """

    # Unpack the input data
    z_coord_list, r_coord_list = profile_data
    width_inj, width_conv, width_throat, width_exit = width_data
    ht_inj, ht_conv, ht_throat, ht_exit = height_data
    beta_inj, beta_conv, beta_throat, beta_exit = angle_data

    # We create the lists of dimensions of the channels by interpolation
    width_list = interp1d(x_chamber_throat_exit, [width_inj, width_conv, width_throat, width_exit], kind='linear')(z_coord_list)
    ht_list = interp1d(x_chamber_throat_exit, [ht_inj, ht_conv, ht_throat, ht_exit], kind='linear')(z_coord_list)
    beta_list = interp1d(x_chamber_throat_exit, [beta_inj, beta_conv, beta_throat, beta_exit], kind='linear')(z_coord_list)

    # The reversed lists go from the end of the nozzle to the injection plate
    z_coord_list_reversed = np.flip(z_coord_list)
    r_coord_list_reversed = np.flip(r_coord_list)
    width_list_reversed = np.flip(width_list)
    ht_list_reversed = np.flip(ht_list)
    beta_list_reversed = np.flip(beta_list)

    # Create and fill the list of alpha angle (angle about the z-axis)
    alpha_list = np.zeros_like(z_coord_list)
    for i, _ in enumerate(z_coord_list):
        alpha_list[i+1] = alpha_list[i] + np.tan(beta_list[i+1]) *\
            abs(z_coord_list_reversed[i]-z_coord_list_reversed[i+1])/r_coord_list_reversed[i+1]

    # We generate the four points of lists corresponding to the four
    # vertices of the cooling channels
    # The central axis of the engine is located along the z-axis
    xA_list = -(r_coord_list_reversed + wall_thickness) * np.sin(np.deg2rad(alpha_list))\
        - (width_list_reversed / 2) * np.cos(np.deg2rad(alpha_list))
    yA_list = (r_coord_list_reversed + wall_thickness) * np.cos(np.deg2rad(alpha_list))\
        - (width_list_reversed / 2) * np.sin(np.deg2rad(alpha_list))
    xB_list = -(r_coord_list_reversed + wall_thickness) * np.sin(np.deg2rad(alpha_list))\
        + (width_list_reversed / 2) * np.cos(np.deg2rad(alpha_list))
    yB_list = (r_coord_list_reversed + wall_thickness) * np.cos(np.deg2rad(alpha_list))\
        + (width_list_reversed / 2) * np.sin(np.deg2rad(alpha_list))

    xC_list = -(r_coord_list_reversed + wall_thickness + ht_list_reversed) * np.sin(np.deg2rad(alpha_list))\
        + (width_list_reversed / 2) * np.cos(np.deg2rad(alpha_list))
    yC_list = (r_coord_list_reversed + wall_thickness + ht_list_reversed) * np.cos(np.deg2rad(alpha_list))\
        + (width_list_reversed / 2) * np.sin(np.deg2rad(alpha_list))
    xD_list = -(r_coord_list_reversed + wall_thickness + ht_list_reversed) * np.sin(np.deg2rad(alpha_list))\
        - (width_list_reversed / 2) * np.cos(np.deg2rad(alpha_list))
    yD_list = (r_coord_list_reversed + wall_thickness + ht_list_reversed) * np.cos(np.deg2rad(alpha_list))\
        - (width_list_reversed / 2) * np.sin(np.deg2rad(alpha_list))

    # We compute the list of channel center points
    x_center_list = -(r_coord_list_reversed + wall_thickness + ht_list_reversed / 2) * np.sin(np.deg2rad(alpha_list))
    y_center_list = (r_coord_list_reversed + wall_thickness + ht_list_reversed / 2) * np.cos(np.deg2rad(alpha_list))

    return

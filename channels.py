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
from mpl_toolkits.mplot3d import Axes3D


def set_axes_equal(ax):
    '''Set 3D plot axes to equal scale (so 1mm is the same on all axes).'''
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    plot_radius = 0.5 * max([x_range, y_range, z_range])
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def project_onto_plane(v, n):
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    return v - np.sum(v * n, axis=-1, keepdims=True) * n


def generate_channels(profile_data, width_data, height_data, angle_data, wall_thickness,
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
    beta_list_reversed = np.flip(beta_list)

    # Create and fill the list of alpha angle (angle about the z-axis)
    alpha_list_reversed = np.zeros_like(z_coord_list_reversed)
    for i, _ in enumerate(z_coord_list_reversed[:-1]):
        alpha_list_reversed[i+1] = np.rad2deg(np.deg2rad(alpha_list_reversed[i]) + np.tan(np.deg2rad(beta_list_reversed[i+1])) *
                                              abs(z_coord_list_reversed[i]-z_coord_list_reversed[i+1])/r_coord_list_reversed[i+1])
    alpha_list = np.flip(alpha_list_reversed)

    # We generate the four points of lists corresponding to the four vertices of the cooling channels
    # The central axis of the engine is located along the z-axis
    xA_list = -(r_coord_list + wall_thickness) * np.sin(np.deg2rad(alpha_list))\
        - (width_list / 2) * np.cos(np.deg2rad(alpha_list))
    yA_list = (r_coord_list + wall_thickness) * np.cos(np.deg2rad(alpha_list))\
        - (width_list / 2) * np.sin(np.deg2rad(alpha_list))
    xB_list = -(r_coord_list + wall_thickness) * np.sin(np.deg2rad(alpha_list))\
        + (width_list / 2) * np.cos(np.deg2rad(alpha_list))
    yB_list = (r_coord_list + wall_thickness) * np.cos(np.deg2rad(alpha_list))\
        + (width_list / 2) * np.sin(np.deg2rad(alpha_list))

    xC_list = -(r_coord_list + wall_thickness + ht_list) * np.sin(np.deg2rad(alpha_list))\
        + (width_list / 2) * np.cos(np.deg2rad(alpha_list))
    yC_list = (r_coord_list + wall_thickness + ht_list) * np.cos(np.deg2rad(alpha_list))\
        + (width_list / 2) * np.sin(np.deg2rad(alpha_list))
    xD_list = -(r_coord_list + wall_thickness + ht_list) * np.sin(np.deg2rad(alpha_list))\
        - (width_list / 2) * np.cos(np.deg2rad(alpha_list))
    yD_list = (r_coord_list + wall_thickness + ht_list) * np.cos(np.deg2rad(alpha_list))\
        - (width_list / 2) * np.sin(np.deg2rad(alpha_list))

    # We compute the list of channel center points
    x_center_list = -(r_coord_list + wall_thickness + ht_list / 2) * np.sin(np.deg2rad(alpha_list))
    y_center_list = (r_coord_list + wall_thickness + ht_list / 2) * np.cos(np.deg2rad(alpha_list))

    # Compute channel inclination (angle with x-y plane)
    dx = np.diff(x_center_list)
    dy = np.diff(y_center_list)
    dz = np.diff(z_coord_list)
    channel_inclination = np.zeros_like(z_coord_list)
    channel_inclination[:-1] = np.rad2deg(np.arccos(abs(dz)/np.sqrt(dx**2 + dy**2 + dz**2)))
    channel_inclination[-1] = channel_inclination[-2]  # repeat last value

    # Compute channel cross-section area along x-y plane
    AB = np.sqrt((xB_list - xA_list)**2 + (yB_list - yA_list)**2)
    BC = np.sqrt((xC_list - xB_list)**2 + (yC_list - yB_list)**2)
    channel_cross_section = AB * BC

    # Compute effective flow cross-sectional area (normal to flow direction)
    effective_channel_cross_section = channel_cross_section * np.abs(np.cos(np.deg2rad(channel_inclination)))

    # Compute effective wetted perimeter and hydraulic diameter
    # Get centerline direction vectors (unit vectors)
    directions = np.stack([dx, dy, dz], axis=1)
    norms = np.linalg.norm(directions, axis=1)
    directions_unit = np.zeros_like(directions)
    directions_unit = directions / norms[:, None]
    # Repeat last direction for last point
    directions_unit = np.vstack([directions_unit, directions_unit[-1]])
    # For each section, get the 3D coordinates of the corners
    points = np.stack([
        np.column_stack([xA_list, yA_list, z_coord_list]),
        np.column_stack([xB_list, yB_list, z_coord_list]),
        np.column_stack([xC_list, yC_list, z_coord_list]),
        np.column_stack([xD_list, yD_list, z_coord_list])
    ], axis=1)  # shape (N, 4, 3)
    # Compute side vectors
    AB = points[:, 1] - points[:, 0]
    BC = points[:, 2] - points[:, 1]
    CD = points[:, 3] - points[:, 2]
    DA = points[:, 0] - points[:, 3]
    sides = [AB, BC, CD, DA]
    # Project each side onto plane normal to flow direction
    perimeter = np.zeros(len(z_coord_list))
    for side in sides:
        proj = project_onto_plane(side, directions_unit)
        perimeter += np.linalg.norm(proj, axis=1)
    effective_wetted_perimeter = perimeter

    # Hydraulic diameter
    hydraulic_diameter = 4 * effective_channel_cross_section / effective_wetted_perimeter

    # Compute fin thickness along the engine (no inclination)
    fin_thickness = ((2 * np.pi * (r_coord_list + wall_thickness)) - (nbc * width_list)) / nbc

    # Compute fin thickness with inclination
    effective_fin_thickness = fin_thickness*np.cos(np.deg2rad(beta_list))

    # Compute total channel length
    channel_total_length = np.sum(np.sqrt(dx**2 + dy**2 + dz**2))
    ######################### PLOTTING #########################

    # 3D plot of channel centerline and surface of revolution
    fig = plt.figure(dpi=figure_dpi)
    ax = fig.add_subplot(111, projection='3d')
    # Plot channel centerline
    ax.plot(x_center_list, y_center_list, z_coord_list, color='blue', label='Channel centerline')
    # Plot surface of revolution
    theta = np.linspace(0, 2 * np.pi, 60)
    Z, Theta = np.meshgrid(z_coord_list_reversed, theta)
    R = np.tile(r_coord_list_reversed, (len(theta), 1))
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)
    ax.plot_surface(X, Y, Z, color='red', alpha=0.3, linewidth=0, antialiased=True)
    ax.set_xlabel('x (radial) [m]')
    ax.set_ylabel('y (radial) [m]')
    ax.set_zlabel('z (engine axis) [m]')
    ax.set_title('3D Channel Centerline and Engine Surface')
    ax.legend()
    set_axes_equal(ax)
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, channel_inclination, label='Channel Inclination (deg)')
    plt.xlabel('z (engine axis) [m]')
    plt.ylabel('Inclination [deg]')
    plt.title('Channel Inclination vs Engine Axis')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, channel_cross_section, label='Initial CS')
    plt.plot(z_coord_list, effective_channel_cross_section, label='Effective CS')
    plt.xlabel('z (engine axis) [m]')
    plt.ylabel('Cross Section [m²]')
    plt.title('Channel Cross Section vs Engine Axis')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, hydraulic_diameter, label='Hydraulic diameter (m)')
    plt.xlabel('z (engine axis) [m]')
    plt.ylabel('Diameter [m]]')
    plt.title('Hydraulic diameter vs Engine Axis')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, fin_thickness, label='Initial fin thickness (m)')
    plt.plot(z_coord_list, effective_fin_thickness, label='Effective fin thickness (m)')
    plt.xlabel('z (engine axis) [m]')
    plt.ylabel('Thickness [m]]')
    plt.title('Fin thickness vs Engine Axis')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    channel_vertices = {"A": np.column_stack((xA_list, yA_list, z_coord_list)),
                        "B": np.column_stack((xB_list, yB_list, z_coord_list)),
                        "C": np.column_stack((xC_list, yC_list, z_coord_list)),
                        "D": np.column_stack((xD_list, yD_list, z_coord_list))}

    channel_centerline = np.column_stack((x_center_list, y_center_list, z_coord_list))

    return channel_vertices, channel_centerline, channel_inclination, width_list, ht_list, \
        effective_channel_cross_section, hydraulic_diameter, \
        effective_fin_thickness, alpha_list, beta_list, channel_total_length

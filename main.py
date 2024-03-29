"""
Created on Fri Nov 27 14:47:27 2020

Original author: Julien S

Refactored and improved by Mehdi D, Paul B, Paul M, Eve X and Antoine R
"""

import time
import csv
import os
import sys

# Calculations
import numpy as np
import cte_tools as t
from main_solver_notqdm import mainsolver

# Data
from Canaux import canaux

# Graphics
from heatequationsolve import carto2D
from volume3d import carto3d, view3d
import matplotlib.pyplot as plt
from tqdm import tqdm  # For progress bars

# GUI

import tkinter as tk
import cte_gui

from plot_cache import reset_cache


class MainProcess:

    def __init__(self):
        pass

    def process(self, entry_dict):

        # print(entry_dict)
        reset_cache("plot_cache")

        print(
            "██████████████████████████████████ START ███████████████████████████████████")
        print(
            "█                                                                          █")
        print(
            "█ Initialisation                                                           █")
        print(
            "█                                                                          █")

        plot_dir = "plot_cache"

        start_time = time.perf_counter()  # Beginning of the timer
        # %% Initial definitions

        # Name of the engine
        engine_name = str(entry_dict["engine_name"])
        # Distance between two points of calculation
        mesh_size = str(entry_dict["mesh_size"])
        # X coordinates of an engine
        x_coords_filename = f"input/{mesh_size}/{engine_name}/x.txt"
        # Y coordinates of an engine
        y_coords_filename = f"input/{mesh_size}/{engine_name}/y.txt"
        # Engine's parameters (found with CEA)
        input_CEA_data = "input/" + entry_dict["input_CEA_data"]

        # Constant input_data_list
        size2 = 7  # Used for the height of the display in 3D view
        limitation = 0.05  # used to build the scales in 3D view
        figure_dpi = 150  # Dots Per Inch (DPI) for all figures (lower=faster)
        # 0=No plots; 1=Important plots; 2=Less important plots: 3=All plots
        plot_detail = int(entry_dict["plot_detail"])
        show_3d_plots = bool(entry_dict["show_3d_plots"])
        show_2D_temperature = bool(entry_dict["show_2D_temperature"])
        do_final_3d_plot = bool(entry_dict["do_final_3d_plot"])
        write_in_csv = bool(entry_dict["write_in_csv"])

        if plot_detail >= 1:
            plot_dir1 = os.path.join(plot_dir, "Basics")
            if not os.path.exists(plot_dir1) and not os.path.isdir(plot_dir1):
                os.mkdir(plot_dir1)

        if plot_detail >= 2:
            plot_dir2 = os.path.join(plot_dir, "Advanced")
            if not os.path.exists(plot_dir2) and not os.path.isdir(plot_dir2):
                os.mkdir(plot_dir2)

        if plot_detail >= 3:
            plot_dir3 = os.path.join(plot_dir, "Global")
            if not os.path.exists(plot_dir3) and not os.path.isdir(plot_dir3):
                os.mkdir(plot_dir3)

        # %% Reading input data
        input_data_reader = csv.reader(open(input_CEA_data, "r"))
        input_data_list = [row[1] for row in input_data_reader]

        # Store CEA output in lists
        # Sound velocity in the chamber
        sound_speed_init = float(input_data_list[0])
        # Sound velocity in the throat
        sound_speed_throat = float(input_data_list[1])
        debit_LOX = float(input_data_list[2])  # LOX debit
        debit_mass_coolant = float(input_data_list[3])  # CH4 debit
        rho_init = float(input_data_list[4])  # Initial density of the gases
        Pc = float(input_data_list[5])  # Pressure in the chamber
        Tc = float(input_data_list[6])  # Combustion temperature
        gamma_c_input = float(input_data_list[7])  # Gamma in the chamber
        gamma_t_input = float(input_data_list[8])  # Gamma in the throat
        gamma_e_input = float(input_data_list[9])  # Gamma at the exit
        molar_mass = float(input_data_list[10])  # Molar mass of the gases
        c_star = float(input_data_list[11])  # Caracteristic velocity
        # Molar fraction of the H2O in the chamber
        xH2O_c_input = float(input_data_list[16])
        # Molar fraction of the H20 in the throat
        xH2O_t_input = float(input_data_list[17])
        # Molar fraction of the H20 at the exit
        xH2O_e_input = float(input_data_list[18])
        # Molar fraction of the CO2 in the chamber
        xCO2_c_input = float(input_data_list[19])
        # Molar fraction of the CO2 in the throat
        xCO2_t_input = float(input_data_list[20])
        # Molar fraction of the CO2 at the exit
        xCO2_e_input = float(input_data_list[21])

        # Store input dimensions in lists
        # Radius of curvature before the throat
        curv_radius_pre_throat = float(input_data_list[12])
        # Radius of curvature after the throat
        curv_radius_after_throat = float(input_data_list[13])
        area_throat = float(input_data_list[14])  # Area at the throat
        diam_throat = float(input_data_list[15])  # Throat diameter

        # Detailed gamma
        gamma_1_input = float(input_data_list[22])
        gamma_2_input = float(input_data_list[23])
        gamma_3_input = float(input_data_list[24])
        gamma_4_input = float(input_data_list[25])
        gamma_5_input = float(input_data_list[26])
        gamma_6_input = float(input_data_list[27])
        gamma_7_input = float(input_data_list[28])
        gamma_8_input = float(input_data_list[29])
        gamma_9_input = float(input_data_list[30])
        gamma_10_input = float(input_data_list[31])
        gamma_11_input = float(input_data_list[32])
        gamma_12_input = float(input_data_list[33])
        gamma_13_input = float(input_data_list[34])
        gamma_14_input = float(input_data_list[35])

        x_1_input = float(input_data_list[36])
        x_2_input = float(input_data_list[37])
        x_3_input = float(input_data_list[38])
        x_4_input = float(input_data_list[39])
        x_5_input = float(input_data_list[40])
        x_6_input = float(input_data_list[41])
        x_7_input = float(input_data_list[42])
        x_8_input = float(input_data_list[43])
        x_9_input = float(input_data_list[44])
        x_10_input = float(input_data_list[45])
        x_11_input = float(input_data_list[46])
        x_12_input = float(input_data_list[47])
        x_13_input = float(input_data_list[48])
        x_14_input = float(input_data_list[49])

        # %% Import of the (X,Y) coordinates of the Viserion
        x_coords_reader = csv.reader(open(x_coords_filename, "r"))
        y_coords_reader = csv.reader(open(y_coords_filename, "r"))

        # Storing the X,Y coordinates in lists
        x_coord_list = [float(row[0]) / 1000 for row in x_coords_reader]
        y_coord_list = [float(row[0]) / 1000 for row in y_coords_reader]
        x_coords_reader = [x * 1000 for x in x_coord_list]
        y_coords_reader = [y * 1000 for y in y_coord_list]

        # Number of points (or the index of the end of the divergent)
        nb_points = len(x_coord_list)

        # Plot of the profile of the engine
        if plot_detail >= 3:
            fig = t.one_plot(x=x_coords_reader, y=y_coords_reader,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Radius of the chamber [mm]',
                             title="Engine profile", show=False,
                             equal_axes=True, ymin=0, ymax=50,
                             xmin=-190, xmax=35)
            fig.savefig(f"{plot_dir3}/profile.png")
            plt.close(fig)

        # Computation and plot of the mesh density of the engine
        if plot_detail >= 3 and show_3d_plots:
            dist_between_pts = [abs(x_coord_list[i] - x_coord_list[i + 1])
                                for i in range(0, len(x_coord_list) - 1)]
            dist_between_pts.append(dist_between_pts[-1])
            colormap = plt.cm.binary
            inv = 1, 1, 1  # 1 means should be reversed
            view3d(inv, x_coord_list, y_coord_list, dist_between_pts,
                   colormap, 'Mesh density (in m)', size2, limitation, plot_dir)

        # %% Computation of the cross-sectional area along the engine
        cross_section_area_list = [np.pi * r ** 2 for r in y_coord_list]

        # Plots of the cross-sectionnal areas
        if plot_detail >= 3:
            fig = t.one_plot(x=x_coords_reader, y=cross_section_area_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Area [m²]',
                             title="Cross-sectional area", show=False,
                             xmin=-190, xmax=35,
                             sci_notation=True, ymin=0)
            fig.savefig(f"{plot_dir3}/Cross section area.png")
            plt.close(fig)

        # %% Adiabatic constant (gamma) parametrization
        print(
            "█ Computing gamma                                                          █")

        i_conv = 0  # Index of the beginning of the convergent
        y1 = 1
        y2 = 1
        while y1 == y2:  # Read y values two per two in order to detect the beginning of the convergent
            y1 = y_coord_list[i_conv]
            i_conv += 1
            y2 = y_coord_list[i_conv]
        i_throat = y_coord_list.index(min(y_coord_list))  # Throat index

        # Gamma in the cylindrical chamber
        # Gamma is constant before the beginning of the convergent
        x_given = [x_coord_list[0],
                   x_coord_list[i_conv],
                   x_1_input,
                   x_coord_list[i_throat],
                   x_2_input,
                   x_3_input,
                   x_4_input,
                   x_5_input,
                   x_6_input,
                   x_7_input,
                   x_8_input,
                   x_9_input,
                   x_10_input,
                   x_11_input,
                   x_12_input,
                   x_13_input,
                   x_14_input,
                   x_coord_list[-1]]
        gamma_given = [gamma_c_input,
                       gamma_c_input,
                       gamma_1_input,
                       gamma_t_input,
                       gamma_2_input,
                       gamma_3_input,
                       gamma_4_input,
                       gamma_5_input,
                       gamma_6_input,
                       gamma_7_input,
                       gamma_8_input,
                       gamma_9_input,
                       gamma_10_input,
                       gamma_11_input,
                       gamma_12_input,
                       gamma_13_input,
                       gamma_14_input,
                       gamma_e_input]
        gamma_list = [x for x in np.interp(x_coord_list, x_given, gamma_given)]

        # Plot of the gamma linearisation
        if plot_detail >= 3:
            fig = t.one_plot(x=x_coords_reader, y=gamma_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Gamma [-]',
                             title="Gamma of hot gases as a function of engine axis",
                             show=False, xmin=-190, xmax=35)
            fig.savefig(f"{plot_dir3}/Gamma of hot gases.png")
            plt.close(fig)

        # %% Mach number computation
        "Computation of gases mach number of the hot gases (and their initial velocity)"

        mach_list = [0.07]

        # Mach number computations along the engine

        # with tqdm(total=nb_points - 1,
        #           desc="█ Computing mach number        ",
        #           unit="|   █", bar_format="{l_bar}{bar}{unit}",
        #           ncols=76) as progressbar:
        #     for i in range(0, nb_points - 1):
        #         mach_gas = t.mach_solv(cross_section_area_list[i], cross_section_area_list[i + 1],
        #                                mach_gas, gamma_list[i])
        #         mach_list.append(mach_gas)
        #         progressbar.update(1)
        print(
            "█ Computing mach number                                                    █")
        for i in range(0, nb_points - 1):
            subsonic = True if i < i_throat else False
            mach_gas = t.mach_solv(cross_section_area_list[i], cross_section_area_list[i_throat],
                                   gamma_list[i], subsonic=subsonic)
            mach_list.append(mach_gas)

        # Plots of the Mach number in the engine (2D/3D)
        if plot_detail >= 1:
            fig = t.one_plot(x=x_coords_reader, y=mach_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Mach [-]',
                             title="Mach number as a function of the engine axis",
                             show=False, xmin=-190, xmax=35)
            fig.savefig(f"{plot_dir3}/Mach number.png")
            plt.close(fig)

        if plot_detail >= 1 and show_3d_plots:
            colormap = plt.cm.Spectral
            inv = 1, 1, 1  # 1 means should be reversed
            print(
                "█ Plotting 3D graph                                                        █")
            view3d(inv, x_coord_list, y_coord_list, mach_list, colormap,
                   'Mach number of hot gases', size2, limitation, plot_dir)

        # %% Static pressure computation
        pressure_list = [Pc]  # (in Pa)

        # with tqdm(total=nb_points - 1,
        #           desc="█ Computing static pressure    ",
        #           unit="|   █", bar_format="{l_bar}{bar}{unit}",
        #           ncols=76) as progressbar:
        #     for i in range(0, nb_points - 1):
        #         pressure = t.pressure_solv(
        #             mach_list[i], mach_list[i + 1], pressure_list[i], gamma_list[i])
        #         pressure_list.append(pressure)
        #         progressbar.update(1)
        print(
            "█ Computing static pressure                                                █")
        for i in range(0, nb_points - 1):
            pressure = t.pressure_solv(
                mach_list[i], mach_list[i + 1], pressure_list[i], gamma_list[i])
            pressure_list.append(pressure)

        # Plot of the static pressure (2D/3D)
        if plot_detail >= 2:
            fig = t.one_plot(x=x_coords_reader, y=pressure_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Pressure [Pa]',
                             title="Global static pressure",
                             show=False, xmin=-190, xmax=35,
                             sci_notation=True, ymin=0)
            fig.savefig(f"{plot_dir2}/Global static pressure.png")
            plt.close(fig)

        if plot_detail >= 2 and show_3d_plots:
            colormap = plt.cm.gist_rainbow_r
            inv = 1, 1, 1  # 1 means should be reversed
            print(
                "█ Plotting 3D graph                                                        █")
            view3d(inv, x_coord_list, y_coord_list, pressure_list,
                   colormap, 'Static pressure (in Pa)', size2, limitation, plot_dir)

        # %% Partial pressure computation and interpolation of the molar fraction
        x_Molfrac = [x_coord_list[0], x_coord_list[i_throat],
                     x_coord_list[-1]]  # Location associated to the molar mass

        # Value of the molar fraction of the H20 after interpolation
        Molfrac_H2O = np.interp(x_coord_list, x_Molfrac, [
            xH2O_c_input, xH2O_t_input, xH2O_e_input])

        # Value of the molar fraction of the CO2 after interpolation
        Molfrac_CO2 = np.interp(x_coord_list, x_Molfrac, [
            xCO2_c_input, xCO2_t_input, xCO2_e_input])

        PH2O_list = [pressure_list[i] * Molfrac_H2O[i]
                     for i in range(0, nb_points)]  # Partial pressure of the H2O
        PCO2_list = [pressure_list[i] * Molfrac_CO2[i]
                     for i in range(0, nb_points)]  # Partial pressure of the CO2

        # Plots of molar fraction and partial pressure
        if plot_detail >= 3:
            fig = t.n_plots(x=x_coords_reader, y_list=[Molfrac_H2O, Molfrac_CO2],
                            xlabel=r'Position x along the engine [mm]',
                            ylabel=r'Molar fraction [-]',
                            y_label_list=["H2O", "CO2"],
                            colors_list=["blue", "orange"],
                            title="Molar fractions of H2O and CO2",
                            show=False, xmin=-190, xmax=35, ymin=0, ymax=1)
            fig.savefig(f"{plot_dir3}/Molar fraction.png")
            plt.close(fig)

            fig = t.n_plots(x=x_coords_reader, y_list=[PH2O_list, PCO2_list],
                            xlabel=r'Position x along the engine [mm]',
                            ylabel=r'Pressure [Pa]',
                            y_label_list=["H2O", "CO2"],
                            colors_list=["blue", "orange"],
                            title="Partial static pressure of H2O and CO2",
                            show=False, xmin=-190, xmax=35, sci_notation=True)
            fig.savefig(f"{plot_dir3}/Partial static pressure.png")
            plt.close(fig)

        # %% Hot gas temperature computation
        hotgas_temp_list = [Tc]
        # with tqdm(total=nb_points - 1,
        #           desc="█ Computing gas static_temperature    ",
        #           unit="|   █", bar_format="{l_bar}{bar}{unit}",
        #           ncols=76) as progressbar:
        #     for i in range(0, nb_points - 1):
        #         static_temperature = t.temperature_hotgas_solv(
        #             mach_list[i], mach_list[i + 1], hotgas_recovery_temp_list[i], gamma_list[i])
        #         hotgas_recovery_temp_list.append(static_temperature)
        #         progressbar.update(1)
        print(
            "█ Computing gas static temperature                                         █")
        for i in range(0, nb_points - 1):
            temperature = t.temperature_hotgas_solv(
                mach_list[i], mach_list[i + 1], hotgas_temp_list[i], gamma_list[i])
            hotgas_temp_list.append(temperature)

        # List of corrected gas temperatures (max diff with original is about 75 K)
        hotgas_temp_list = [t.tempcorrige_pempie(hotgas_temp_list[i], gamma_list[i], mach_list[i]) for i in
                            range(0, nb_points)]

        # Plots of the temperature in the engine (2D/3D)
        if plot_detail >= 2:
            fig = t.one_plot(x=x_coords_reader, y=hotgas_temp_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Temperature [K]',
                             title="Gas static temperature",
                             show=False, xmin=-190, xmax=35)
            fig.savefig(f"{plot_dir2}/Gas temperature.png")
            plt.close(fig)

        if plot_detail >= 2 and show_3d_plots:
            colormap = plt.cm.coolwarm
            inv = 1, 1, 1  # 1 means should be reversed
            print(
                "█ Plotting 3D graph                                                        █")
            view3d(inv, x_coord_list, y_coord_list, hotgas_temp_list, colormap, 'Temperature of the gases (in K)',
                   size2,
                   limitation, plot_dir)

        # %% Dimensions
        print(
            "█ Computing channel geometry                                               █")
        print(
            "█                                                                          █")

        nbc = int(entry_dict["nbc"])  # Number of channels
        # Position of the manifold from the throat (in m)
        manifold_pos = float(entry_dict["manifold_pos"])

        # Widths
        # Width of the channel in at the injection plate (in m)
        lrg_inj = float(entry_dict["lrg_inj"])
        # Width of the channel at the end of the cylindrical chamber (in m)
        lrg_conv = float(entry_dict["lrg_conv"])
        # Width of the channel in the throat (in m)
        lrg_col = float(entry_dict["lrg_col"])
        # Width of the channel at the manifold (in m)
        lrg_tore = float(entry_dict["lrg_tore"])

        # Heights
        # Height of the channel at the injection plate (in m)
        ht_inj = float(entry_dict["ht_inj"])
        # Height of the channel at the end of the cylindrical chamber (in m)
        ht_conv = float(entry_dict["ht_conv"])
        # Height of the channel in the throat (in m)
        ht_col = float(entry_dict["ht_col"])
        # Height of the channel at the manifold (in m)
        ht_tore = float(entry_dict["ht_tore"])

        # Thickness
        # Thickness of the wall at the chamber (in m)
        e_conv = float(entry_dict["e_conv"])
        # Thickness of the wall at the throat (in m)
        e_col = float(entry_dict["e_col"])
        # Thickness of the wall at the manifold (in m)
        e_tore = float(entry_dict["e_tore"])

        n1 = float(entry_dict["n1"])  # Width convergent
        n2 = float(entry_dict["n2"])  # Width divergent
        n3 = float(entry_dict["n3"])  # Height convergent
        n4 = float(entry_dict["n4"])  # Height divergent
        n5 = float(entry_dict["n5"])  # Thickness convergent
        n6 = float(entry_dict["n6"])  # Thickness divergent

        # %% Material selection
        material_name = entry_dict["material_name"]

        # %% Properties of the coolant
        fluid = entry_dict["fluid"]
        # Density of the CH4 (kg/m^3)
        density_cool_init = float(entry_dict["density_cool_init"])
        # Initial temperature of the coolant (K)
        Temp_cool_init = float(entry_dict["Temp_cool_init"])
        # Total volumic flow rate of the coolant (m^3/s)
        debit_volumique_total_cool = debit_mass_coolant / density_cool_init
        # Pressure of the coolant at inlet (Pa)
        Pressure_cool_init = float(entry_dict["Pressure_cool_init"])
        roughness = float(entry_dict["roughness"])  # Roughness (m)

        # %% Computation of channel geometry

        # Pack the data in tuples
        profile = (x_coord_list, y_coord_list)
        widths = (lrg_inj, lrg_conv, lrg_col, lrg_tore)
        heights = (ht_inj, ht_conv, ht_col, ht_tore)
        thicknesses = (e_conv, e_col, e_tore)
        coeffs = (n1, n2, n3, n4, n5, n6)

        # Compute dimensions
        xcanaux, ycanaux, larg_canal, larg_ailette_list, ht_canal, wall_thickness, area_channel, nb_points_channel, \
        y_coord_avec_canaux \
            = canaux(profile, widths, heights, thicknesses, coeffs, manifold_pos, debit_volumique_total_cool, nbc,
                     plot_detail, write_in_csv, figure_dpi, plot_dir)

        # Write the dimensions of the channels in a CSV file
        file_name = "output/channelvalue.csv"
        with open(file_name, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(("Engine x", "Engine y", "y coolant wall", "Channel width", "Rib width",
                             "Channel height", "Chamber wall thickness", "Channel area"))
            for i in range(0, nb_points_channel):
                writer.writerow((xcanaux[i], y_coord_avec_canaux[i], ycanaux[i], larg_canal[i], larg_ailette_list[i],
                                 ht_canal[i], wall_thickness[i], area_channel[i]))

        end_init_time = time.perf_counter()  # End of the initialisation timer
        # Initialisation elapsed time (in s)
        time_elapsed = f"{round(end_init_time - start_time, 2)}"
        if len(time_elapsed) <= 3:
            time_elapsed_i = f"   {time_elapsed} s"
        elif len(time_elapsed) == 4:
            time_elapsed_i = f"  {time_elapsed} s"
        elif len(time_elapsed) == 5:
            time_elapsed_i = f" {time_elapsed} s"
        else:
            time_elapsed_i = f"{time_elapsed} s"
        start_main_time = time.perf_counter()  # Start of the main solution timer

        # %% Prepare the lists before main computation

        wall_thickness.reverse()
        xcanaux.reverse()
        larg_canal.reverse()
        area_channel.reverse()
        ht_canal.reverse()
        ycanaux.reverse()
        y_coord_avec_canaux.reverse()
        # We reverse the lists in order to calculate from the manifold to the injection

        # Save the data for exporting, before altering the original lists
        hotgas_temperature_saved = hotgas_temp_list[:]
        aire_saved = cross_section_area_list[:]
        mach_list_saved = mach_list[:]
        gamma_saved = gamma_list[:]
        PH2O_list_saved = PH2O_list[:]
        PCO2_list_saved = PCO2_list[:]

        # Remove the data points before the manifold
        hotgas_temp_list = hotgas_temp_list[:nb_points_channel]
        cross_section_area_list = cross_section_area_list[:nb_points_channel]
        mach_list = mach_list[:nb_points_channel]
        gamma_list = gamma_list[:nb_points_channel]
        PH2O_list = PH2O_list[:nb_points_channel]
        PCO2_list = PCO2_list[:nb_points_channel]

        gamma_list.reverse()
        mach_list.reverse()
        cross_section_area_list.reverse()
        hotgas_temp_list.reverse()
        PH2O_list.reverse()
        PCO2_list.reverse()

        # %% Main computation

        data_hotgas = (hotgas_temp_list, molar_mass, gamma_list,
                       Pc, c_star, PH2O_list, PCO2_list)
        data_coolant = (Temp_cool_init, Pressure_cool_init,
                        fluid, debit_mass_coolant)
        data_channel = (xcanaux, ycanaux, larg_canal, larg_ailette_list, ht_canal,
                        wall_thickness, area_channel, nb_points_channel)
        data_chamber = (y_coord_avec_canaux, nbc, diam_throat, curv_radius_pre_throat, area_throat,
                        roughness, cross_section_area_list, mach_list, material_name)

        # Call the main solving loop
        hlcor_list, hlcor_list_2, hotgas_visc_list, hotgas_cp_list, hotgas_cond_list, \
        hotgas_prandtl_list, hg_list, hotwall_temp_list, coldwall_temp_list, flux_list, \
        sigma_list, coolant_reynolds_list, tempcoolant_list, visccoolant_list, \
        condcoolant_list, cpcoolant_list, densitycoolant_list, velocitycoolant_list, \
        pcoolant_list, wallcond_list, sound_speed_coolant_list, hlnormal_list, \
        qRad_list, q_list_CO2, q_list_H2O \
            = mainsolver(data_hotgas, data_coolant, data_channel, data_chamber)

        end_m = time.perf_counter()  # End of the main solution timer
        # Main computation elapsed time (in s)
        time_elapsed = f"{round(end_m - start_main_time, 2)}"
        if len(time_elapsed) <= 3:
            time_elapsed_m = f"   {time_elapsed} s"
        elif len(time_elapsed) == 4:
            time_elapsed_m = f"  {time_elapsed} s"
        elif len(time_elapsed) == 5:
            time_elapsed_m = f" {time_elapsed} s"
        else:
            time_elapsed_m = f"{time_elapsed} s"

        # %% Display of the 1D analysis results
        print(
            "█                                                                          █")

        if plot_detail >= 1:
            start_d1 = time.perf_counter()  # Start of the display of 1D timer
            print(
                "█ Displaying results                                                       █")
            print(
                "█                                                                          █")
            fig = t.n_plots(x=x_coords_reader,
                            y_list=[hlcor_list_2, hlcor_list, hlnormal_list],
                            xlabel=r'Position x along the engine [mm]',
                            ylabel=r'Convection coeff [W/m² K]',
                            y_label_list=["Hl corrected (Luka Denies)",
                                          "Hl corrected (Julien)",
                                          "Hl normal"],
                            colors_list=["blue", "green", "cyan"],
                            title="Cold-side convection coefficient",
                            show=False, xmin=-190, xmax=35, reverse=True,
                            ymin=0)
            fig.savefig(f"{plot_dir1}/Convection coeff.png")
            plt.close(fig)

            fig = t.n_plots(x=x_coords_reader,
                            y_list=[coldwall_temp_list, hotwall_temp_list],
                            xlabel=r'Position x along the engine [mm]',
                            ylabel=r'Temperature [K]',
                            y_label_list=["Twl", "Twg"],
                            colors_list=["blue", "red"],
                            title="Wall temperatures",
                            show=False, xmin=-190, xmax=35,
                            reverse=True, ymin=273)
            fig.savefig(f"{plot_dir1}/Wall temperature.png")
            plt.close(fig)

            tempcoolant_list.pop()
            fig = t.one_plot(x=x_coords_reader, y=tempcoolant_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Temperature [K]',
                             title="Coolant temperature",
                             show=False, xmin=-190, xmax=35, reverse=True)
            fig.savefig(f"{plot_dir1}/Coolant temperature.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=flux_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Total heat flux density [W/m²]',
                             title="Heat flux density",
                             show=False, xmin=-190, xmax=35,
                             sci_notation=True, reverse=True, ymin=0)
            fig.savefig(f"{plot_dir1}/Heat flux.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=velocitycoolant_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Velocity [m/s]',
                             title="Coolant velocity",
                             show=False, xmin=-190, xmax=35,
                             reverse=True)
            fig.savefig(f"{plot_dir1}/Velocity.png")
            plt.close(fig)

            pcoolant_list.pop()
            fig = t.one_plot(x=x_coords_reader, y=pcoolant_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Pressure [Pa]',
                             title="Pressure drop in the cooling channels",
                             show=False, xmin=-190, xmax=35,
                             sci_notation=True, reverse=True, ymin=0)
            fig.savefig(f"{plot_dir1}/Pressure drop.png")
            plt.close(fig)

        if plot_detail >= 2:
            fig = t.one_plot(x=x_coords_reader, y=wallcond_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Conductivity [W/m K]',
                             title="Conductivity of the wall",
                             show=False, xmin=-190, xmax=35, reverse=True)
            fig.savefig(f"{plot_dir2}/Wall conductivity.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=hg_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Convection coefficient [W/m² K]',
                             title="Gas-side convection coefficient",
                             show=False, xmin=-190, xmax=35, reverse=True,
                             ymin=0)
            fig.savefig(f"{plot_dir2}/Convection coefficient.png")
            plt.close(fig)

            densitycoolant_list.pop()
            fig = t.one_plot(x=x_coords_reader, y=densitycoolant_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Density [$\frac{kg}{m^3}$]',
                             title="Coolant density",
                             show=False, xmin=-190, xmax=35,
                             reverse=True)
            fig.savefig(f"{plot_dir2}/Coolant Volumic mass.png")
            plt.close(fig)

            fig = t.n_plots(x=x_coords_reader,
                            y_list=[q_list_CO2, q_list_H2O, qRad_list],
                            xlabel=r'Position x along the engine [mm]',
                            ylabel=r'Heat flux density [W/m²]',
                            y_label_list=["CO2", "H2O", "Total"],
                            colors_list=["r", "b", "k"],
                            title="Radiative heat flux density",
                            show=False, xmin=-190, xmax=35, ymin=0,
                            sci_notation=True, reverse=True)
            fig.savefig(f"{plot_dir2}/Radiative heat flux.png")
            plt.close(fig)

        if plot_detail >= 3:
            fig = t.one_plot(x=x_coords_reader, y=coolant_reynolds_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Reynolds [-]',
                             title="Coolant Reynolds number",
                             show=False, xmin=-190, xmax=35, reverse=True)
            fig.savefig(f"{plot_dir3}/Reynolds number.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=hotgas_visc_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Dynamic viscosity [Pa s]',
                             title="Coolant viscosity",
                             show=False, xmin=-190, xmax=35,
                             sci_notation=True, reverse=True)
            fig.savefig(f"{plot_dir3}/Gas viscosity.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=hotgas_cp_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Cp [J/kg K]',
                             title="Gas specific heat",
                             show=False, xmin=-190, xmax=35, reverse=True)
            fig.savefig(f"{plot_dir3}/Gas Cp.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=hotgas_cond_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Conductivity [W/m K]',
                             title="Gas conductivity",
                             show=False, xmin=-190, xmax=35,
                             sci_notation=True, reverse=True)
            fig.savefig(f"{plot_dir3}/Gas conductivity.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=hotgas_prandtl_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Prandtl [-]',
                             title="Gas Prandtl number",
                             show=False, xmin=-190, xmax=35, reverse=True)
            fig.savefig(f"{plot_dir3}/Gas Prandtl.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=sigma_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Sigma [-]',
                             title="Bartz correction factor sigma",
                             show=False, xmin=-190, xmax=35, reverse=True)
            fig.savefig(f"{plot_dir3}/Sigma.png")
            plt.close(fig)

            condcoolant_list.pop()
            fig = t.one_plot(x=x_coords_reader, y=condcoolant_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Conductivity [W/m K]',
                             title="Coolant conductivity",
                             show=False, xmin=-190, xmax=35,
                             sci_notation=True, reverse=True)
            fig.savefig(f"{plot_dir3}/Coolant conductivity.png")
            plt.close(fig)

            cpcoolant_list.pop()
            fig = t.one_plot(x=x_coords_reader, y=cpcoolant_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Specific heat [J/kg K]',
                             title="Coolant specific heat",
                             show=False, xmin=-190, xmax=35, reverse=True)
            fig.savefig(f"{plot_dir3}/Coolant Cp.png")
            plt.close(fig)

            visccoolant_list.pop()
            fig = t.one_plot(x=x_coords_reader, y=visccoolant_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Dynamic viscosity [Pa s]',
                             title="Coolant viscosity",
                             show=False, xmin=-190, xmax=35,
                             sci_notation=True, reverse=True)
            fig.savefig(f"{plot_dir3}/Coolant viscosity.png")
            plt.close(fig)

            fig = t.one_plot(x=x_coords_reader, y=sound_speed_coolant_list,
                             xlabel=r'Position x along the engine [mm]',
                             ylabel=r'Sound speed [m/s]',
                             title="Coolant speed of sound",
                             show=False, xmin=-190, xmax=35, reverse=True)
            fig.savefig(f"{plot_dir3}/Sound velocity.png")
            plt.close(fig)

        if plot_detail >= 1 and show_3d_plots:
            colormap = plt.cm.plasma
            inv = 0, 0, 0

            view3d(inv, xcanaux, ycanaux, flux_list, colormap,
                   "Heat flux (in MW m2)", size2, limitation, plot_dir)

            colormap = plt.cm.coolwarm
            inv = 0, 0, 0
            view3d(inv, xcanaux, ycanaux, tempcoolant_list, colormap,
                   "Temperature of the coolant (in K)", size2, limitation, plot_dir)

        if plot_detail >= 2 and show_3d_plots:
            colormap = plt.cm.magma
            inv = 0, 0, 0  # 1 means should be reversed
            view3d(inv, xcanaux, ycanaux, coldwall_temp_list, colormap,
                   "Wall temperature on the gas side (in K)", size2, limitation, plot_dir)

        if plot_detail >= 1:
            end_d1 = time.perf_counter()  # End of the display of 1D timer
            # 1D display elapsed time (in s)
            time_elapsed = f"{round(end_d1 - start_d1, 2)}"
            if len(time_elapsed) <= 3:
                time_elapsed_d1 = f"   {time_elapsed} s"
            elif len(time_elapsed) == 4:
                time_elapsed_d1 = f"  {time_elapsed} s"
            elif len(time_elapsed) == 5:
                time_elapsed_d1 = f" {time_elapsed} s"
            else:
                time_elapsed_d1 = f"{time_elapsed} s"

        # %% Flux computation in 2D and 3D
        """2D flux computation"""
        larg_ailette_list.reverse()

        if show_2D_temperature:
            start_d2 = time.perf_counter()  # Start of the display of 2D timer
            # At the beginning of the chamber
            print(
                "█ Results at the beginning of the chamber :                                █")
            dx = 0.00004  # *3.5
            location = " at the beginning of the chamber"
            carto2D(larg_ailette_list[-1] + larg_canal[-1], larg_canal[-1], e_conv, ht_canal[-1], dx, hg_list[-1],
                    wallcond_list[-1], hotgas_temp_list[-1], hlcor_list[-1], tempcoolant_list[-1], 5, True, 1, location,
                    False, plot_dir)

            # At the throat
            print(
                "█ Results at the throat :                                                  █")
            pos_col = ycanaux.index(min(ycanaux))
            dx = 0.000025  # *3.5
            location = " at the throat"
            carto2D(larg_ailette_list[pos_col] + larg_canal[pos_col], larg_canal[pos_col], e_col, ht_canal[pos_col],
                    dx, hg_list[pos_col], wallcond_list[pos_col], hotgas_temp_list[pos_col], hlcor_list[pos_col],
                    tempcoolant_list[pos_col], 15, True, 2, location, False, plot_dir)
            # At the end of the divergent
            print(
                "█ Results at the manifold :                                                █")
            dx = 0.00004
            location = " at the manifold"
            carto2D(larg_ailette_list[0] + larg_canal[0], larg_canal[0], e_tore, ht_canal[0], dx, hg_list[0],
                    wallcond_list[0], hotgas_temp_list[0], hlcor_list[0], tempcoolant_list[0], 5, True, 1, location,
                    False, plot_dir)

            end_d2 = time.perf_counter()  # End of the display of 2D timer
            # 2D display elapsed time (in s)
            time_elapsed = f"{round(end_d2 - start_d2, 2)}"
            if len(time_elapsed) <= 3:
                time_elapsed_d2 = f"   {time_elapsed} s"
            elif len(time_elapsed) == 4:
                time_elapsed_d2 = f"  {time_elapsed} s"
            elif len(time_elapsed) == 5:
                time_elapsed_d2 = f" {time_elapsed} s"
            else:
                time_elapsed_d2 = f"{time_elapsed} s"

        "Computation for 3D graph"
        if do_final_3d_plot:
            start_d3 = time.perf_counter()
            temperature_slice_list = []
            lim1 = 0
            lim2 = 650
            dx = 0.0001

            # Compute a (low-resolution) 2D slice for each point in the engine
            """with tqdm(total=nb_points_channel,
                    desc="█ 3D graph computation         ",
                    unit="|   █", bar_format="{l_bar}{bar}{unit}",
                    ncols=76) as progressbar:
                for i in range(0, nb_points_channel):
                    temperature_slice = carto2D(larg_ailette_list[i] + larg_canal[i], larg_canal[i], wall_thickness[i],
                                                ht_canal[i], dx, hg_list[i], wallcond_list[i], 
                                                hotgas_recovery_temp_list[i],
                                                hlnormal_list[i], tempcoolant_list[i], 3, False, 1, "", True)
                    temperature_slice_list.append(temperature_slice)
                    progressbar.update(1)"""

            print("█ 3D graph computation")
            for i in range(0, nb_points_channel):
                temperature_slice = carto2D(larg_ailette_list[i] + larg_canal[i], larg_canal[i], wall_thickness[i],
                                            ht_canal[i], dx, hg_list[i], wallcond_list[i], hotgas_temp_list[i],
                                            hlnormal_list[i], tempcoolant_list[i], 3, False, 1, "", True, plot_dir)
                temperature_slice_list.append(temperature_slice)

            # Stack all these slices in a final 3D plot
            carto3d([0, 0, 0], xcanaux, ycanaux, temperature_slice_list, plt.cm.Spectral_r,
                    '3D view of wall temperatures (in K)', nbc, limitation, plot_dir)
            print(
                "█                                                                          █")
            # End of the 3D display timer
            end_d3 = time.perf_counter()
            # 3D display elapsed time (in s)
            time_elapsed = f"{round(end_d3 - start_d3, 2)}"
            if len(time_elapsed) <= 3:
                time_elapsed_d3 = f"   {time_elapsed} s"
            elif len(time_elapsed) == 4:
                time_elapsed_d3 = f"  {time_elapsed} s"
            elif len(time_elapsed) == 5:
                time_elapsed_d3 = f" {time_elapsed} s"
            else:
                time_elapsed_d3 = f"{time_elapsed} s"
        start_e = time.perf_counter()  # Start of the end timer

        # %% Reversion of the lists

        cross_section_area_list.reverse()
        gamma_list.reverse()
        mach_list.reverse()
        hotgas_temp_list.reverse()
        xcanaux.reverse()
        ycanaux.reverse()
        larg_canal.reverse()
        ht_canal.reverse()
        area_channel.reverse()
        hotgas_visc_list.reverse()
        hotgas_cp_list.reverse()
        hotgas_cond_list.reverse()
        hg_list.reverse()
        hotgas_prandtl_list.reverse()
        sigma_list.reverse()
        coldwall_temp_list.reverse()
        hotwall_temp_list.reverse()
        flux_list.reverse()
        tempcoolant_list.reverse()
        velocitycoolant_list.reverse()
        coolant_reynolds_list.reverse()
        hlnormal_list.reverse()
        densitycoolant_list.reverse()
        visccoolant_list.reverse()
        condcoolant_list.reverse()
        cpcoolant_list.reverse()
        pcoolant_list.reverse()
        PH2O_list.reverse()
        PCO2_list.reverse()
        y_coord_avec_canaux.reverse()

        # %% Preparation of the lists for CAD modelisation
        "Changing the coordinates of the height of the channels (otherwise it is geometrically wrong)"

        angles = [0]
        newxhtre = [xcanaux[0]]
        newyhtre = [ycanaux[0] + ht_canal[0]]
        for i in range(1, nb_points_channel):
            if i == (nb_points_channel - 1):
                angle = angles[i - 1]
                angles.append(angle)
            else:
                vect1 = (xcanaux[i] - xcanaux[i - 1]) / (
                        (((ycanaux[i] - ycanaux[i - 1]) ** 2) + ((xcanaux[i] - xcanaux[i - 1]) ** 2)) ** 0.5)
                vect2 = (xcanaux[i + 1] - xcanaux[i]) / (
                        (((ycanaux[i + 1] - ycanaux[i]) ** 2) + ((xcanaux[i + 1] - xcanaux[i]) ** 2)) ** 0.5)
                angle1 = np.rad2deg(np.arccos(vect1))
                angle2 = np.rad2deg(np.arccos(vect2))
                angle = angle2
                angles.append(angle)
            newx = xcanaux[i] + ht_canal[i] * np.sin(np.deg2rad(angles[i]))
            newy = ycanaux[i] + ht_canal[i] * np.cos(np.deg2rad(angles[i]))
            newxhtre.append(newx)
            newyhtre.append(newy)

        # Checking the height of channels
        verification = []
        print(
            "█ Checking and computing channel height                                    █")
        for i in range(0, nb_points_channel):
            verifhtre = (((newxhtre[i] - xcanaux[i]) ** 2) +
                         ((newyhtre[i] - ycanaux[i]) ** 2)) ** 0.5
            verification.append(verifhtre)

        if plot_detail >= 3:
            plt.figure(dpi=figure_dpi)
            plt.plot(newxhtre, newyhtre, color='blue', label='New height')
            plt.plot(xcanaux, ycanaux, color='chocolate',
                     label='Former height')
            plt.title(
                "Geometrical aspect of the channel (height as a function of the engine axis)")
            plt.axis("equal")
            plt.legend(loc='upper left')
            plt.savefig(f"{plot_dir3}/Geometrical aspect.png")
            plt.close()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, verification)
            plt.title("Checking the height of the generated channels")
            plt.savefig(f"{plot_dir3}/Height checking.png")
            plt.close()

        # %% Writing the results of the study in a CSV file

        if write_in_csv:
            print(
                "█ Writing results in .csv files                                            █")
            valuexport = open("output/valuexport.csv", "w", newline="")
            geometry1 = open("output/geometry1.csv", "w", newline="")
            geometry2 = open("output/geometry2.csv", "w", newline="")
            valuexport_writer = csv.writer(valuexport)
            geometry1_writer = csv.writer(geometry1)
            geometry2_writer = csv.writer(geometry2)
            valuexport_writer.writerow(
                ("Engine x axix", "Engine diameter", "Area of gas engine", "Gas gamma",
                 "Mach number", "Gas pressure", "Total pressure", "Gas static temperature",
                 "Channels x axis", "Engine + chamber wall diameter", "Channels width",
                 "Channels height", "Channels area", "Gas viscosity", "Cp gas",
                 "Gas conductivity", "Prandtl gaz", "Coeff Hg", "Sigma", " Twg ", " Twl ",
                 "Heat flux", "Tl", "Reynolds CH4", "Coeff Hl", "Rho coolant",
                 "Viscosity CH4", "Conductivity CH4", "Cp CH4", "Coolant velocity",
                 "Coolant pressure", "Wall conductivity", "x real height", "y real height"))
            geometry1_writer.writerow(("x real height", "y real height"))
            geometry2_writer.writerow(
                ("Engine + chamber wall radius", "x real height"))

            for i in range(0, nb_points):
                if i < nb_points_channel:
                    geometry1_writer.writerow(
                        (newxhtre[i] * (-1000), newyhtre[i] * 1000))
                    geometry2_writer.writerow(
                        (ycanaux[i] * 1000, newxhtre[i] * (-1000)))
                    valuexport_writer.writerow((x_coord_list[i], y_coord_list[i], aire_saved[i], gamma_saved[i],
                                                mach_list_saved[i], pressure_list[i],
                                                hotgas_temperature_saved[i], xcanaux[i], ycanaux[i],
                                                larg_canal[i], ht_canal[i], area_channel[i], hotgas_visc_list[i],
                                                hotgas_cp_list[i], hotgas_cond_list[i], hotgas_prandtl_list[i],
                                                hg_list[i], sigma_list[i], coldwall_temp_list[i],
                                                hotwall_temp_list[i], flux_list[i], tempcoolant_list[i],
                                                coolant_reynolds_list[i], hlnormal_list[i], densitycoolant_list[i],
                                                visccoolant_list[i],
                                                condcoolant_list[i], cpcoolant_list[i], velocitycoolant_list[i],
                                                pcoolant_list[i],
                                                wallcond_list[i], newxhtre[i], newyhtre[i]))
                else:
                    valuexport_writer.writerow(
                        (x_coord_list[i], y_coord_list[i], aire_saved[i], gamma_saved[i],
                         mach_list_saved[i], pressure_list[i],
                         hotgas_temperature_saved[i], ' ', ' ', ' ', ' ', ' ', ' ',
                         ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                         ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))

            valuexport.close()
            geometry1.close()
            geometry2.close()

            longc = len(xcanaux)

            # write in a csv for CATIA macro (work in progress)
            if write_in_csv:
                "Writing the results of the study in a CSV file"
                file_name = "output/channel_macro_catia.csv"
                file = open(file_name, "w", newline="")
                writer = csv.writer(file)
                writer.writerow(["StartCurve"])
                for i in range(0, longc, 3):
                    writer.writerow(
                        (1000 * xcanaux[i], 1000 * (ycanaux[i]), 1000 * (larg_canal[i] / 2)))
                writer.writerow(["EndCurve"])
                writer.writerow(["StartCurve"])
                for i in range(0, longc, 3):
                    writer.writerow(
                        (1000 * xcanaux[i], 1000 * (ycanaux[i]), 1000 * (-larg_canal[i] / 2)))
                writer.writerow(["EndCurve"])
                writer.writerow(["StartCurve"])
                for i in range(0, longc, 3):
                    writer.writerow(
                        (1000 * xcanaux[i], 1000 * (ycanaux[i] + ht_canal[i]), 1000 * (larg_canal[i] / 2)))
                writer.writerow(["EndCurve"])
                writer.writerow(["StartCurve"])
                for i in range(0, longc, 3):
                    writer.writerow(
                        (1000 * xcanaux[i], 1000 * (ycanaux[i] + ht_canal[i]), 1000 * (- larg_canal[i] / 2)))
                writer.writerow(["EndCurve"])

                # connection between two adjacent edges (exit)
                writer.writerow(["StartCurve"])
                writer.writerow((1000 * xcanaux[0], 1000 * ycanaux[0], 500 * larg_canal[0]))
                writer.writerow((1000 * xcanaux[0], 1000 * ycanaux[0], -500 * larg_canal[0]))
                writer.writerow(["EndCurve"])

                writer.writerow(["StartCurve"])
                writer.writerow((1000 * xcanaux[0], 1000 * ycanaux[0], 500 * larg_canal[0]))
                writer.writerow((1000 * xcanaux[0], 1000 * (ycanaux[0] + ht_canal[0]), 500 * larg_canal[0]))
                writer.writerow(["EndCurve"])

                writer.writerow(["StartCurve"])
                writer.writerow((1000 * xcanaux[0], 1000 * (ycanaux[0] + ht_canal[0]), 500 * larg_canal[0]))
                writer.writerow((1000 * xcanaux[0], 1000 * (ycanaux[0] + ht_canal[0]), -500 * larg_canal[0]))
                writer.writerow(["EndCurve"])

                writer.writerow(["StartCurve"])
                writer.writerow((1000 * xcanaux[0], 1000 * (ycanaux[0] + ht_canal[0]), -500 * larg_canal[0]))
                writer.writerow((1000 * xcanaux[0], 1000 * ycanaux[0], -500 * larg_canal[0]))
                writer.writerow(["EndCurve"])

                # connection between two adjacent edges (injection)
                writer.writerow(["StartCurve"])
                writer.writerow((1000 * xcanaux[-1], 1000 * ycanaux[-1], 500 * larg_canal[-1]))
                writer.writerow((1000 * xcanaux[-1], 1000 * ycanaux[-1], -500 * larg_canal[-1]))
                writer.writerow(["EndCurve"])

                writer.writerow(["StartCurve"])
                writer.writerow((1000 * xcanaux[-1], 1000 * ycanaux[-1], 500 * larg_canal[-1]))
                writer.writerow((1000 * xcanaux[-1], 1000 * (ycanaux[-1] + ht_canal[-1]), 500 * larg_canal[-1]))
                writer.writerow(["EndCurve"])

                writer.writerow(["StartCurve"])
                writer.writerow((1000 * xcanaux[-1], 1000 * (ycanaux[-1] + ht_canal[-1]), 500 * larg_canal[-1]))
                writer.writerow((1000 * xcanaux[-1], 1000 * (ycanaux[-1] + ht_canal[-1]), -500 * larg_canal[-1]))
                writer.writerow(["EndCurve"])

                writer.writerow(["StartCurve"])
                writer.writerow((1000 * xcanaux[-1], 1000 * (ycanaux[-1] + ht_canal[-1]), -500 * larg_canal[-1]))
                writer.writerow((1000 * xcanaux[-1], 1000 * ycanaux[-1], -500 * larg_canal[-1]))
                writer.writerow(["EndCurve"])

                writer.writerow(["End"])
                file.close()

        # %% Execution time display

        end_t = time.perf_counter()  # End of the total timer
        # End elapsed time (in s)
        time_elapsed = f"{round(end_t - start_e, 2)}"
        if len(time_elapsed) <= 3:
            time_elapsed_e = f"   {time_elapsed} s"
        elif len(time_elapsed) == 4:
            time_elapsed_e = f"  {time_elapsed} s"
        elif len(time_elapsed) == 5:
            time_elapsed_e = f" {time_elapsed} s"
        else:
            time_elapsed_e = f"{time_elapsed} s"

        time_elapsed = f"{round(end_t - start_time, 2)}"  # Total elapsed time
        if len(time_elapsed) <= 3:
            time_elapsed_t = f"   {time_elapsed} s"
        elif len(time_elapsed) == 4:
            time_elapsed_t = f"  {time_elapsed} s"
        elif len(time_elapsed) == 5:
            time_elapsed_t = f" {time_elapsed} s"
        else:
            time_elapsed_t = f"{time_elapsed} s"

        print(
            "█                                                                          █")
        print(
            "█__________________________________________________________________________█")
        print(
            "█                                                                          █")
        print(
            f"█ Execution time for the initialisation       : {time_elapsed_i}                   █")
        print(
            "█                                                                          █")
        print(
            f"█ Execution time for the main computation     : {time_elapsed_m}                   █")

        if plot_detail >= 1:
            print(
                "█                                                                          █")
            print(
                f"█ Execution time for the display of 1D        : {time_elapsed_d1}                   █")

        if show_2D_temperature:
            print(
                "█                                                                          █")
            print(
                f"█ Execution time for the display of 2D        : {time_elapsed_d2}                   █")

        if do_final_3d_plot:
            print(
                "█                                                                          █")
            print(
                f"█ Execution time for the display of 3D        : {time_elapsed_d3}                   █")

        print(
            "█                                                                          █")
        print(
            f"█ Execution time for the end of the program   : {time_elapsed_e}                   █")
        print(
            "█                                                                          █")
        print(
            f"█ Total execution time                        : {time_elapsed_t}                   █")
        print(
            "█                                                                          █")
        print(
            "███████████████████████████████████ END ████████████████████████████████████")


plot_dir = "plot_cache"
reset_cache(plot_dir)

gui = cte_gui.MainGUI(MainProcess)
gui.title("CTE")

"""input_class = cte_gui.InputsWin()
entry_dict = input_class.entry_dict"""

# runprocess = cte_gui.Run(MainProcess)


print("██████████████████████████ Cool The Engine V 2.0.0 █████████████████████████")
print("█                                                                          █")
print("█                  Innovative Propulsion Laboratory - IPL                  █")
print("█__________________________________________________________________________█")

while True:
    try:
        gui.update()
        gui.mainloop()
        break

    except tk.TclError:
        break
# print("ok")
reset_cache("plot_cache")
exit()

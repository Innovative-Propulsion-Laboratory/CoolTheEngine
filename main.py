"""
Created on Fri Nov 27 14:47:27 2020

Original author: Julien S

Refactored and improved by Mehdi D, Paul B, Paul M, Eve X and Antoine R
"""

import time
import csv

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


class MainProcess:
    def __init__(self):
        pass

    def process(self):
        print(
            "██████████████████████████████████ START ███████████████████████████████████")
        print(
            "█                                                                          █")
        print(
            "█ Initialisation                                                           █")
        print(
            "█                                                                          █")

        start_time = time.perf_counter()  # Beginning of the timer
        # %% Initial definitions

        mesh_size = 0.25  # Distance between two points of calculation
        # X coordinates of the Viserion
        x_coords_filename = f"input/{mesh_size}/x.txt"
        # Y coordinates of the Viserion
        y_coords_filename = f"input/{mesh_size}/y.txt"
        # Viserion's parameters (found with CEA)
        input_CEA_data = "input/Viserion_2023.txt"

        # Constant input_data_list
        size2 = 16  # Used for the height of the display in 3D view
        limitation = 0.05  # used to build the scales in 3D view
        figure_dpi = 150  # Dots Per Inch (DPI) for all figures (lower=faster)
        plot_detail = 0  # 0=No plots; 1=Important plots; 2=Less important plots: 3=All plots
        show_3d_plots = False
        show_2D_temperature = False
        do_final_3d_plot = False
        write_in_csv = False

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

        # %% Import of the (X,Y) coordinates of the Viserion
        x_coords_reader = csv.reader(open(x_coords_filename, "r"))
        y_coords_reader = csv.reader(open(y_coords_filename, "r"))

        # Storing the X,Y coordinates in lists
        x_coord_list = [float(row[0]) / 1000 for row in x_coords_reader]
        y_coord_list = [float(row[0]) / 1000 for row in y_coords_reader]
        # Number of points (or the index of the end of the divergent)
        nb_points = len(x_coord_list)

        # Plot of the profile of the engine
        if plot_detail >= 3:
            plt.figure(dpi=figure_dpi)
            plt.plot(x_coord_list, y_coord_list, color='black')
            plt.title(
                'Profile of the Viserion (left : chamber and right : divergent)', color='black')
            plt.show()

        # Computation and plot of the mesh density of the engine
        if plot_detail >= 3 and show_3d_plots:
            dist_between_pts = [abs(x_coord_list[i] - x_coord_list[i + 1])
                                for i in range(0, len(x_coord_list) - 1)]
            dist_between_pts.append(dist_between_pts[-1])
            colormap = plt.cm.binary
            inv = 1, 1, 1  # 1 means should be reversed
            view3d(inv, x_coord_list, y_coord_list, dist_between_pts,
                   colormap, 'Mesh density (in m)', size2, limitation)

        # %% Computation of the cross-sectional area along the engine
        cross_section_area_list = [np.pi * r ** 2 for r in y_coord_list]

        # Plots of the cross-sectionnal areas
        if plot_detail >= 3:
            plt.figure(dpi=figure_dpi)
            plt.plot(x_coord_list, cross_section_area_list, color='black')
            plt.title(
                "Cross section area of engine (in m²) as a function of engine axis")
            plt.show()

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

        # Gamma in the cylindrical chamber
        # Gamma is constant before the beginning of the convergent
        gamma_list = [gamma_c_input for i in range(0, i_conv)]

        # Gamma in the convergent
        i_throat = y_coord_list.index(min(y_coord_list))  # Throat index
        gamma_convergent = gamma_c_input
        for m in range(-1, i_throat - i_conv - 1):
            # Linear interpolation between beginning and end of convergent:
            # (yi+1)=((y2-y1)/(x2-x1))*abs((xi+1)-(xi))
            gamma_convergent += ((gamma_t_input - gamma_c_input) / (x_coord_list[i_throat] - x_coord_list[i_conv])) * abs(
                x_coord_list[i_conv + 1 + m] - x_coord_list[i_conv + m])
            gamma_list.append(gamma_convergent)

        # Gamma in the divergent nozzle
        gamma_divergent = gamma_t_input
        # Linear interpolation between beginning and end of divergent
        for q in range(-1, nb_points - i_throat - 1):
            gamma_divergent += ((gamma_e_input - gamma_t_input) / (x_coord_list[-1] - x_coord_list[i_throat])) * abs(
                x_coord_list[i_throat + 1 + q] - x_coord_list[i_throat + q])
            gamma_list.append(gamma_divergent)

        # Plot of the gamma linearisation
        if plot_detail >= 3:
            plt.figure(dpi=figure_dpi)
            plt.plot(x_coord_list, gamma_list, color='gold')
            plt.title("Gamma of hot gases as a function of engine axis")
            plt.show()

        # %% Mach number computation
        "Computation of gases mach number of the hot gases (and their initial velocity)"

        v_init_gas = (debit_LOX + debit_mass_coolant) / (rho_init *
                                                         cross_section_area_list[0])  # Initial velocity of the gases
        mach_init_gas = v_init_gas / sound_speed_init  # Initial mach number
        mach_gas = mach_init_gas
        mach_list = [mach_init_gas]

        # Mach number computations along the engine
        """with tqdm(total=nb_points - 1,
                desc="█ Computing mach number        ",
                unit="|   █", bar_format="{l_bar}{bar}{unit}",
                ncols=76) as progressbar:
            for i in range(0, nb_points - 1):
                mach_gas = t.mach_solv(cross_section_area_list[i], cross_section_area_list[i + 1],
                                    mach_gas, gamma_list[i])
                mach_list.append(mach_gas)
                progressbar.update(1)"""
        print(
            "█ Computing mach number                                                    █")
        for i in range(0, nb_points - 1):
            mach_gas = t.mach_solv(cross_section_area_list[i], cross_section_area_list[i + 1],
                                   mach_gas, gamma_list[i])
            mach_list.append(mach_gas)

        # Plots of the Mach number in the engine (2D/3D)
        if plot_detail >= 1:
            plt.figure(dpi=figure_dpi)
            plt.plot(x_coord_list, mach_list, color='gold')
            plt.title("Mach number as a function of the engine axis")
            plt.show()

        if plot_detail >= 1 and show_3d_plots:
            colormap = plt.cm.Spectral
            inv = 1, 1, 1  # 1 means should be reversed
            print(
                "█ Plotting 3D graph                                                        █")
            view3d(inv, x_coord_list, y_coord_list, mach_list, colormap,
                   'Mach number of hot gases', size2, limitation)

        # %% Static pressure computation
        pressure_list = [Pc]  # (in Pa)

        """with tqdm(total=nb_points - 1,
                desc="█ Computing static pressure    ",
                unit="|   █", bar_format="{l_bar}{bar}{unit}",
                ncols=76) as progressbar:
            for i in range(0, nb_points - 1):
                pressure = t.pressure_solv(
                    mach_list[i], mach_list[i + 1], pressure_list[i], gamma_list[i])
                pressure_list.append(pressure)
                progressbar.update(1)"""
        print(
            "█ Computing static pressure                                                █")
        for i in range(0, nb_points - 1):
            pressure = t.pressure_solv(
                mach_list[i], mach_list[i + 1], pressure_list[i], gamma_list[i])
            pressure_list.append(pressure)

        # Plot of the static pressure (2D/3D)
        if plot_detail >= 2:
            plt.figure(dpi=figure_dpi)
            plt.plot(x_coord_list, pressure_list, color='gold')
            plt.title(
                "Global static pressure (in Pa) as a function of the engine axis")
            plt.show()

        if plot_detail >= 2 and show_3d_plots:
            colormap = plt.cm.gist_rainbow_r
            inv = 1, 1, 1  # 1 means should be reversed
            print(
                "█ Plotting 3D graph                                                        █")
            view3d(inv, x_coord_list, y_coord_list, pressure_list,
                   colormap, 'Static pressure (in Pa)', size2, limitation)

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
            plt.figure(dpi=figure_dpi)
            plt.plot(x_coord_list, Molfrac_H2O, color='blue', label='H20')
            plt.plot(x_coord_list, Molfrac_CO2, color='orange', label='C02')
            plt.title("Molar fraction of as a function of the engine axis")
            plt.legend(loc='center left')
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(x_coord_list, PH2O_list, color='blue', label='H20')
            plt.plot(x_coord_list, PCO2_list, color='orange', label='C02')
            plt.title(
                "Partial static pressure (in Pa) of as a function of the engine axis")
            plt.legend(loc='center left')
            plt.show()

        # %% Hot gas temperature computation
        hotgas_temp_list = [Tc]
        """with tqdm(total=nb_points - 1,
                desc="█ Computing gas temperature    ",
                unit="|   █", bar_format="{l_bar}{bar}{unit}",
                ncols=76) as progressbar:
            for i in range(0, nb_points - 1):
                temperature = t.temperature_hotgas_solv(
                    mach_list[i], mach_list[i + 1], hotgas_temp_list[i], gamma_list[i])
                hotgas_temp_list.append(temperature)
                progressbar.update(1)"""
        print(
            "█ Computing gas temperature                                                █")
        for i in range(0, nb_points - 1):
            temperature = t.temperature_hotgas_solv(
                mach_list[i], mach_list[i + 1], hotgas_temp_list[i], gamma_list[i])
            hotgas_temp_list.append(temperature)

        # List of corrected gas temperatures (max diff with original is about 75 K)
        hotgas_temp_list = [t.tempcorrige(hotgas_temp_list[i], gamma_list[i], mach_list[i]) for i in
                            range(0, nb_points)]

        # Plots of the temperature in the engine (2D/3D)
        if plot_detail >= 2:
            plt.figure(dpi=figure_dpi)
            plt.plot(x_coord_list, hotgas_temp_list, color='gold')
            plt.title("Gas temperature (in K) as a function of the engine axis")
            plt.show()

        if plot_detail >= 2 and show_3d_plots:
            colormap = plt.cm.coolwarm
            inv = 1, 1, 1  # 1 means should be reversed
            print(
                "█ Plotting 3D graph                                                        █")
            view3d(inv, x_coord_list, y_coord_list, hotgas_temp_list, colormap, 'Temperature of the gases (in K)', size2,
                   limitation)

        # %% Dimensions
        print(
            "█ Computing channel geometric                                              █")
        print(
            "█                                                                          █")

        nbc = 40  # Number of channels
        manifold_pos = 0.104  # Position of the manifold from the throat (in m)

        # Widths
        # Width of the channel in at the injection plate (in m)
        lrg_inj = 0.0045
        # Width of the channel at the end of the cylindrical chamber (in m)
        lrg_conv = 0.0025
        lrg_col = 0.0015  # Width of the channel in the throat (in m)
        lrg_tore = 0.002  # Width of the channel at the manifold (in m)

        # Heights
        ht_inj = 0.002  # Height of the channel at the injection plate (in m)
        # Height of the channel at the end of the cylindrical chamber (in m)
        ht_conv = 0.002
        ht_col = 0.0015  # Height of the channel in the throat (in m)
        ht_tore = 0.002  # Height of the channel at the manifold (in m)

        # Thickness
        e_conv = 0.001  # Thickness of the wall at the chamber (in m)
        e_col = 0.001  # Thickness of the wall at the throat (in m)
        e_tore = 0.001  # Thickness of the wall at the manifold (in m)

        n1 = 1  # Width convergent
        n2 = 1  # Width divergent
        n3 = 1  # Height convergent
        n4 = 1  # Height divergent
        n5 = 1  # Thickness convergent
        n6 = 1  # Thickness divergent

        # %% Material selection
        material = 1
        if material == 0:
            material_name = "pure copper"
        elif material == 1:
            material_name = "cucrzr"
        elif material == 2:
            material_name = "inconel"

        # %% Properties of the coolant
        fluid = "Methane"
        density_cool_init = 425  # Density of the CH4 (kg/m^3)
        Temp_cool_init = 110  # Initial temperature of the coolant (K)
        # Total volumic flow rate of the coolant (m^3/s)
        debit_volumique_total_cool = debit_mass_coolant / density_cool_init
        Pressure_cool_init = 7000000  # Pressure of the coolant at inlet (Pa)
        roughness = 15e-6  # Roughness (m)

        # %% Computation of channel geometry

        # Pack the data in tuples
        profile = (x_coord_list, y_coord_list)
        widths = (lrg_inj, lrg_conv, lrg_col, lrg_tore)
        heights = (ht_inj, ht_conv, ht_col, ht_tore)
        thicknesses = (e_conv, e_col, e_tore)
        coeffs = (n1, n2, n3, n4, n5, n6)

        # Compute dimensions
        xcanaux, ycanaux, larg_canal, larg_ailette_list, ht_canal, wall_thickness, area_channel, nb_points_channel, y_coord_avec_canaux \
            = canaux(profile, widths, heights, thicknesses, coeffs, manifold_pos, debit_volumique_total_cool, nbc, plot_detail,
                     write_in_csv, figure_dpi)

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
                "█ Display of results                                                       █")
            print(
                "█                                                                          █")
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, hlcor_list_2, color='blue',
                     label='Hl corrected (Luka Denies)')
            plt.plot(xcanaux, hlcor_list, color='green',
                     label='Hl corrected (Julien)')
            plt.plot(xcanaux, hlnormal_list, color='cyan', label='Hl')
            plt.title("Convection coeff as a function of the engine axis")
            plt.legend(loc='upper left')
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, coldwall_temp_list, color='blue', label='Twl')
            plt.plot(xcanaux, hotwall_temp_list, color='red', label='Twg')
            plt.title('Wall temperature (in K) as a function of engine axis')
            plt.legend(loc='lower left')
            plt.show()

            tempcoolant_list.pop()
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, tempcoolant_list, color='blue')
            plt.title('Coolant temperature (in K) as a function of engine axis')
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, flux_list, color='red')
            plt.title('Heat flux (in W) as a function of engine axis')
            plt.show()

            mach_03 = [x * 0.3 for x in sound_speed_coolant_list]
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, velocitycoolant_list,
                     color='blue', label='Coolant')
            plt.plot(xcanaux, mach_03, color='orange', label='Mach 0.3 limit')
            plt.title(
                'Velocity (in m/s) of the coolant as a function of engine axis')
            plt.legend(loc='upper left')
            plt.show()

            pcoolant_list.pop()
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, pcoolant_list, color='orange')
            plt.title('Pressure drop in the cooling channels')
            plt.show()

        if plot_detail >= 2:
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, wallcond_list, color='orangered')
            plt.title('Conductivity of the wall as a function of engine axis')
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, hg_list, color='orangered')
            plt.title('Convection coefficient Hg as a function of engine axis')
            plt.show()

            plt.figure(dpi=figure_dpi)
            densitycoolant_list.pop()
            plt.plot(xcanaux, densitycoolant_list, color='blue')
            plt.title('Volumic mass of the coolant as a function of engine axis')
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, q_list_CO2, color='r', label='CO2')
            plt.plot(xcanaux, q_list_H2O, color='b', label='H2O')
            plt.plot(xcanaux, qRad_list, color='g', label='total')
            plt.title('Radiative heat flux(W/m2)')
            plt.legend()
            plt.show()

        if plot_detail >= 3:
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, coolant_reynolds_list, color='blue')
            plt.title(
                "Reynolds number of the coolant as a function of the engine axis")
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, hotgas_visc_list, color='orangered')
            plt.title("Gas viscosity as a function of the engine axis")
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, hotgas_cp_list, color='orangered')
            plt.title("Gas Cp as a function of the engine axis")
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, hotgas_cond_list, color='orangered')
            plt.title("Gas conductivity as a function of engine axis")
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, hotgas_prandtl_list, color='orangered')
            plt.title("Gas Prandtl number as a function of engine axis")
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, sigma_list, color='orangered')
            plt.title("Sigma as a function of the engine axis")
            plt.show()

            condcoolant_list.pop()
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, condcoolant_list, color='blue')
            plt.title('Conductivity of the coolant as a function of engine axis')
            plt.show()

            cpcoolant_list.pop()
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, cpcoolant_list, color='blue')
            plt.title('Cp of the coolant as a function of engine axis')
            plt.show()

            visccoolant_list.pop()
            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, visccoolant_list, color='blue')
            plt.title('Viscosity of the coolant as a function of engine axis')
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, sound_speed_coolant_list, color='pink')
            plt.title(
                'Sound velocity of the coolant (in m/s) as a function of engine axis')
            plt.show()

        if plot_detail >= 1 and show_3d_plots:
            colormap = plt.cm.plasma
            inv = 0, 0, 0
            view3d(inv, xcanaux, ycanaux, flux_list, colormap,
                   "Heat flux (in MW/m²)", size2, limitation)

            colormap = plt.cm.coolwarm
            inv = 0, 0, 0
            view3d(inv, xcanaux, ycanaux, tempcoolant_list, colormap,
                   "Temperature of the coolant (in K)", size2, limitation)

        if plot_detail >= 2 and show_3d_plots:
            colormap = plt.cm.magma
            inv = 0, 0, 0  # 1 means should be reversed
            view3d(inv, xcanaux, ycanaux, coldwall_temp_list, colormap, "Wall temperature on the gas side (in K)", size2,
                   limitation)

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
                    wallcond_list[-1], hotgas_temp_list[-1], hlcor_list[-1], tempcoolant_list[-1], 5, True, 1, location, False)

            # At the throat
            print(
                "█ Results at the throat :                                                  █")
            pos_col = ycanaux.index(min(ycanaux))
            dx = 0.000025  # *3.5
            location = " at the throat"
            carto2D(larg_ailette_list[pos_col] + larg_canal[pos_col], larg_canal[pos_col], e_col, ht_canal[pos_col],
                    dx, hg_list[pos_col], wallcond_list[pos_col], hotgas_temp_list[pos_col], hlcor_list[pos_col],
                    tempcoolant_list[pos_col], 15, True, 2, location, False)
            # At the end of the divergent
            print(
                "█ Results at the manifold :                                                █")
            dx = 0.00004
            location = " at the manifold"
            carto2D(larg_ailette_list[0] + larg_canal[0], larg_canal[0], e_tore, ht_canal[0], dx, hg_list[0],
                    wallcond_list[0], hotgas_temp_list[0], hlcor_list[0], tempcoolant_list[0], 5, True, 1, location, False)

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
                                                ht_canal[i], dx, hg_list[i], wallcond_list[i], hotgas_temp_list[i],
                                                hlnormal_list[i], tempcoolant_list[i], 3, False, 1, "", True)
                    temperature_slice_list.append(temperature_slice)
                    progressbar.update(1)"""

            print("█ 3D graph computation")
            for i in range(0, nb_points_channel):
                temperature_slice = carto2D(larg_ailette_list[i] + larg_canal[i], larg_canal[i], wall_thickness[i],
                                            ht_canal[i], dx, hg_list[i], wallcond_list[i], hotgas_temp_list[i],
                                            hlnormal_list[i], tempcoolant_list[i], 3, False, 1, "", True)
                temperature_slice_list.append(temperature_slice)

            # Stack all these slices in a final 3D plot
            carto3d([0, 0, 0], xcanaux, ycanaux, temperature_slice_list, plt.cm.Spectral_r,
                    '3D view of wall temperatures (in K)', nbc, limitation)
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
            plt.show()

            plt.figure(dpi=figure_dpi)
            plt.plot(xcanaux, verification)
            plt.title("Checking the height of the generated channels")
            plt.show()

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
                 "Mach number", "Gas pressure", "Total pressure", "Gas temperature",
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


gui = cte_gui.MainGUI(MainProcess)
gui.title("CTE")

runprocess = cte_gui.Run(MainProcess)

print("██████████████████████████ Cool The Engine V 2.0.0 █████████████████████████")
print("█                                                                          █")
print("█                  Innovative Propulsion Laboratory - IPL                  █")
print("█__________________________________________________________________________█")


while True:
    try:
        gui.update()
        gui.mainloop()
        exit()
    except tk.TclError:
        break
exit()

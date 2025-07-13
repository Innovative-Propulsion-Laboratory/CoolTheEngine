"""
Created on Sat Jul 12 2025

Original author: Paul M
"""

from scipy.interpolate import interp1d, PchipInterpolator
import csv
import json
from itertools import product
import pandas as pd

# Calculations
import numpy as np
import cte_tools as t
from solver import solver

# Data
from channels import generate_channels
import cea

# Graphics
from plotter import plotter


class TaskManager():
    def __init__(self, json_path: str):
        with open(json_path, "r") as f:
            self.settings = json.load(f)

        # Check if any parameter is a list, then run bulk mode
        if any(isinstance(param, list) for param in list(self.settings.values())):
            self.bulk_run()

        # Otherwise, run the single mode (no parameter is a list)
        else:
            self.single_run()

    def bulk_run(self):
        # Transform into singletons if not a list, and store all parameters in a list of lists
        value_lists = [v if isinstance(v, list) else [v] for v in (self.settings.values())]

        nb_configs = np.prod([len(lst) for lst in value_lists])
        if nb_configs > 5000:
            raise ValueError(f"Too many configurations ({nb_configs}>5000)")

        print(f"Solving {nb_configs} unique configurations.")
        # Generate every combination
        combinations = list(product(*value_lists))

        # Store every configuration in a dataframe
        configurations_df = pd.DataFrame(combinations, columns=(self.settings.keys()))
        configurations_df["cte_instance"] = pd.Series(dtype="object")
        configurations_df["max_wall_temp"] = pd.Series(dtype="float")
        configurations_df["avg_wall_temp"] = pd.Series(dtype="float")
        configurations_df["max_stress"] = pd.Series(dtype="float")
        configurations_df["coolant_pressure_drop"] = pd.Series(dtype="float")
        configurations_df["coolant_temp_increase"] = pd.Series(dtype="float")
        for i, row in configurations_df.iterrows():
            # Create and store the instance of CoolTheEngine
            cte = CoolTheEngine(dict(row))
            configurations_df.loc[i, "cte_instance"] = cte

            # Compute and unpack results
            try:
                (max_wall_temp, avg_wall_temp, max_stress,
                    coolant_pressure_drop, coolant_temp_increase) = cte.compute_all()
            except ValueError as err:
                print(f"Failed computation for config {i}.\t Error: {err}")
                (max_wall_temp, avg_wall_temp, max_stress,
                    coolant_pressure_drop, coolant_temp_increase) = -1, -1, -1, -1, -1

            # Assign results to respective columns
            configurations_df.loc[i, "max_wall_temp"] = max_wall_temp
            configurations_df.loc[i, "avg_wall_temp"] = avg_wall_temp
            configurations_df.loc[i, "max_stress"] = max_stress
            configurations_df.loc[i, "coolant_pressure_drop"] = coolant_pressure_drop
            configurations_df.loc[i, "coolant_temp_increase"] = coolant_temp_increase

        print(configurations_df)
        configurations_df.to_csv("output/bulk_results.csv")

    def single_run(self):
        cte = CoolTheEngine(self.settings)
        cte.compute_all(show_1D=True, show_2D=True, save_plot=True, write_in_csv=True)


class CoolTheEngine():
    def __init__(self, params: dict):
        # Engine parameters
        self.chamber_pressure = float(params["chamber_pressure"])
        self.coolant_inlet_pressure = float(params["P_coolant"])
        self.coolant_inlet_temp = float(params["T_coolant"])
        self.ox_mfr = float(params["ox_mfr"])
        self.fuel_mfr = float(params["fuel_mfr"])
        self.coolant_mfr = float(params["coolant_mfr"])
        self.ox_name = params["oxidizer_name"]
        self.fuel_name = params["fuel_name"]
        self.coolant_name = params["coolant_name"]
        self.use_TEOS_PDMS = params["use_TEOS_PDMS"]

        # Channel parameters
        self.wall_material = params["wall_material"]
        self.channel_roughness = float(params["channel_roughness"])
        self.nb_channels = int(params["nb_channels"])
        self.wall_thickness = float(params["wall_thickness"])

        # Widths
        self.width_inj = float(params["channel_widths_inj"])
        self.width_conv = float(params["channel_widths_conv"])
        self.width_throat = float(params["channel_widths_throat"])
        self.width_exit = float(params["channel_widths_exit"])

        # Heights
        self.ht_inj = float(params["channel_heights_inj"])
        self.ht_conv = float(params["channel_heights_conv"])
        self.ht_throat = float(params["channel_heights_throat"])
        self.ht_exit = float(params["channel_heights_exit"])

        # Angles
        self.beta_inj = float(params["channel_angles_inj"])
        self.beta_conv = float(params["channel_angles_conv"])
        self.beta_throat = float(params["channel_angles_throat"])
        self.beta_exit = float(params["channel_angles_exit"])

        # General options
        self.nb_points = int(params["nb_points"])
        self.figure_dpi = int(params["figure_dpi"])

        # Initial definitions
        self.contour_file = "input/engine_contour.csv"  # Engine contour

    def compute_all(self, show_1D=False, show_2D=False, save_plot=False,
                    write_in_csv=False, mach_list=None):
        # Reading input data
        contour_data = np.genfromtxt(self.contour_file, delimiter=",", skip_header=1)
        z_coord_list = contour_data[0:, 0]/1000
        r_coord_list = contour_data[0:, 1]/1000
        nb_points_raw = len(z_coord_list)  # Number of points

        # Reduce the number of points using scipy 1D interpolation
        # Create new x values evenly spaced between the min and max of the original x_coord_list
        x_new = np.linspace(z_coord_list[0], z_coord_list[-1], self.nb_points)
        # Interpolate y values at the new x positions
        interp_func = interp1d(z_coord_list, r_coord_list, kind='linear')
        y_new = interp_func(x_new)
        # Replace the original lists with the interpolated ones
        z_coord_list = x_new
        r_coord_list = y_new

        throat_radius = np.min(r_coord_list)  # Radius of the throat (in m)
        exit_radius = r_coord_list[-1]  # Radius of the exit (in m)
        throat_diam = 2 * throat_radius  # Diameter of the throat (in m)

        throat_area = np.pi * throat_radius**2  # Area of the throat (in m²)
        exit_area = np.pi * exit_radius**2  # Area of the exit (in m²)

        expansion_ratio = exit_area / throat_area  # Expansion ratio (Ae/A*)
        throat_curv_radius = 1.5 * throat_radius  # Curvature radius before the throat (in m)

        # Find the index of the throat
        i_throat = np.argmin(np.abs(r_coord_list))
        # Find the index of the beginning of the convergent section
        i_convergent = 0
        for i in range(1, len(r_coord_list)):
            if r_coord_list[i] < r_coord_list[i-1]:
                i_convergent = i
                break

        Cstar, Tc, MolWt = cea.compute_Cstar_Tc_MolWt(self.chamber_pressure, self.ox_mfr/self.fuel_mfr,
                                                      self.ox_name, self.fuel_name, expansion_ratio)

        hotgas_mu_chamber, hotgas_cp_chamber, hotgas_lambda_chamber, hotgas_pr_chamber, \
            hotgas_mu_throat, hotgas_cp_throat, hotgas_lambda_throat, hotgas_pr_throat, \
            hotgas_mu_exit, hotgas_cp_exit, hotgas_lambda_exit, hotgas_pr_exit = cea.get_hotgas_properties(self.chamber_pressure, self.ox_mfr/self.fuel_mfr,
                                                                                                           self.ox_name, self.fuel_name, expansion_ratio)

        x_chamber_throat_exit = [z_coord_list[0], z_coord_list[i_convergent],
                                 z_coord_list[i_throat], z_coord_list[-1]]

        hotgas_visc_list = PchipInterpolator(x_chamber_throat_exit,
                                             [hotgas_mu_chamber, hotgas_mu_chamber,
                                              hotgas_mu_throat, hotgas_mu_exit])(z_coord_list)
        hotgas_cp_list = PchipInterpolator(x_chamber_throat_exit,
                                           [hotgas_cp_chamber, hotgas_cp_chamber,
                                            hotgas_cp_throat, hotgas_cp_exit])(z_coord_list)
        hotgas_cond_list = PchipInterpolator(x_chamber_throat_exit,
                                             [hotgas_lambda_chamber, hotgas_lambda_chamber,
                                              hotgas_lambda_throat, hotgas_lambda_exit])(z_coord_list)
        hotgas_pr_list = PchipInterpolator(x_chamber_throat_exit,
                                           [hotgas_pr_chamber, hotgas_pr_chamber,
                                            hotgas_pr_throat, hotgas_pr_exit])(z_coord_list)

        # Computation of the cross-sectional area along the engine
        cross_section_area_list = np.pi*r_coord_list**2

        gamma_list = cea.compute_gamma(self.chamber_pressure, self.ox_mfr/self.fuel_mfr,
                                       self.ox_name, self.fuel_name,
                                       cross_section_area_list/throat_area)

        # Computation of mach number of the hot gases
        mach_list = t.mach_list_from_area_ratios(cross_section_area_list/throat_area,
                                                 gamma_list[0], i_throat)

        # Static pressure computation
        static_pressure_list = np.zeros_like(z_coord_list)  # (in Pa)
        for i in range(0, self.nb_points):
            static_pressure_list[i] = t.pressure_solv(mach_list[i], gamma_list[i],
                                                      self.chamber_pressure)

        # Partial pressure computation and interpolation of the molar fraction
        molFracH2O_chamber, molFracH2O_throat, molFracH2O_exit = cea.compute_H2O_molar_fractions(self.chamber_pressure, self.ox_mfr/self.fuel_mfr,
                                                                                                 self.ox_name, self.fuel_name, expansion_ratio)
        molFracCO2_chamber, molFracCO2_throat, molFracCO2_exit = cea.compute_CO2_molar_fractions(self.chamber_pressure, self.ox_mfr/self.fuel_mfr,
                                                                                                 self.ox_name, self.fuel_name, expansion_ratio)

        # Linear interpolation of molar fractions for H2O and CO2
        molFracH2O = PchipInterpolator(x_chamber_throat_exit, [molFracH2O_chamber,
                                                               molFracH2O_chamber,
                                                               molFracH2O_throat,
                                                               molFracH2O_exit])(z_coord_list)
        molFracCO2 = PchipInterpolator(x_chamber_throat_exit, [molFracCO2_chamber,
                                                               molFracCO2_chamber,
                                                               molFracCO2_throat,
                                                               molFracCO2_exit])(z_coord_list)

        # Partial pressure of H2O and CO2
        P_H2O_list = np.array([static_pressure_list[i] * molFracH2O[i] for i in range(0, self.nb_points)])
        P_CO2_list = np.array([static_pressure_list[i] * molFracCO2[i] for i in range(0, self.nb_points)])

        # Hot gas temperature computation
        hotgas_static_temp_list = np.zeros_like(z_coord_list)  # List of static hot gas temperatures
        for i in range(0, self.nb_points):
            hotgas_static_temp_list[i] = t.temperature_hotgas_solv(mach_list[i], gamma_list[i], Tc)

        # We assume that the total temperature is constant
        hotgas_total_temp_list = Tc*np.ones_like(z_coord_list)  # List of total hot gas temperatures

        # Computation of adiabatic wall temperature (recovery temperature)
        hotgas_recovery_temp_list = np.array([t.get_recovery_temperature(hotgas_total_temp_list[i], gamma_list[i], mach_list[i], hotgas_pr_list[i]) for i in
                                              range(0, self.nb_points)])

        # Pack the data in tuples
        profile = (z_coord_list, r_coord_list)
        widths = (self.width_inj, self.width_conv, self.width_throat, self.width_exit)
        heights = (self.ht_inj, self.ht_conv, self.ht_throat, self.ht_exit)
        angles = (self.beta_inj, self.beta_conv, self.beta_throat, self.beta_exit)

        # Generate the cooling channels
        channel_vertices, channel_centerline, channel_inclination, channel_ar_list, \
            channel_width_list, channel_height_list, \
            initial_channel_cross_section, effective_channel_cross_section, \
            hydraulic_diameter, initial_fin_thickness, effective_fin_thickness_list, \
            alpha_list, beta_list, channel_total_length = generate_channels(profile, widths, heights,
                                                                            angles, self.wall_thickness,
                                                                            self.nb_channels, x_chamber_throat_exit)

        # Main computation
        data_hotgas = (hotgas_recovery_temp_list, hotgas_static_temp_list, hotgas_visc_list, hotgas_pr_list,
                       hotgas_cp_list, hotgas_cond_list, MolWt, gamma_list, self.chamber_pressure,
                       Cstar, P_H2O_list, P_CO2_list)
        data_coolant = (self.coolant_inlet_temp, self.coolant_inlet_pressure,
                        self.coolant_name, self.coolant_mfr, self.use_TEOS_PDMS)
        data_channel = (self.nb_channels, channel_width_list,
                        channel_height_list, effective_fin_thickness_list,
                        self.wall_thickness, hydraulic_diameter,
                        effective_channel_cross_section,
                        channel_centerline, beta_list)
        data_chamber = (z_coord_list, r_coord_list, throat_diam, throat_curv_radius, throat_area,
                        self.channel_roughness, cross_section_area_list, mach_list, self.wall_material)

        # Call the main solving loop
        hl_corrected_list, h_tp_list, hg_list, \
            hotwall_temp_list, coldwall_temp_list, q_tot_list, sigma_list, \
            coolant_reynolds_list, coolant_temp_list, coolant_visc_list, \
            coolant_cond_list, coolant_cp_list, coolant_density_list, \
            coolant_velocity_list, coolant_pressure_list, coolant_Tsat_list, wall_cond_list, \
            hl_normal_list, hl_corrected_list, q_rad_list, q_rad_list_CO2, q_rad_list_H2O, \
            CHF_Meyer_list, CHF_Tong_list = solver(data_hotgas, data_coolant, data_channel, data_chamber)

        hoop_stress_list, thermal_stress_list, max_wall_stress_list = t.compute_1D_wall_stress(self.wall_material, self.wall_thickness, z_coord_list,
                                                                                               r_coord_list, static_pressure_list,
                                                                                               coolant_pressure_list, hotwall_temp_list,
                                                                                               coldwall_temp_list)

        max_wall_temp = np.max(hotwall_temp_list)
        avg_wall_temp = np.average(hotwall_temp_list)
        max_stress = np.max(max_wall_stress_list)
        coolant_pressure_drop = coolant_pressure_list[-1] - coolant_pressure_list[0]
        coolant_temp_increase = coolant_temp_list[0] - coolant_temp_list[-1]

        # PLot the results
        if show_1D or show_2D:
            # Display of the 1D analysis results
            parameters_plotter = (show_1D, show_2D, self.figure_dpi, save_plot)

            # Store the data in a big tuple to send to the plotter
            data_plotter = (  # Engine geometry
                z_coord_list * 1000,
                r_coord_list * 1000,
                z_coord_list,
                r_coord_list,
                cross_section_area_list,

                # Channel geometry
                initial_fin_thickness,
                effective_fin_thickness_list,
                channel_ar_list,
                channel_width_list,
                channel_height_list,
                hydraulic_diameter,
                initial_channel_cross_section,
                effective_channel_cross_section,
                self.nb_channels,
                alpha_list,
                beta_list,
                channel_centerline[:, 0],
                channel_centerline[:, 1],

                # Hotgas properties
                gamma_list,
                mach_list,
                static_pressure_list,
                hotgas_total_temp_list,
                hotgas_recovery_temp_list,
                hotgas_static_temp_list,
                hotgas_visc_list,
                hotgas_cp_list,
                hotgas_cond_list,
                hotgas_pr_list,

                # Convective and radiative heat flux
                hg_list,
                sigma_list,
                hl_normal_list,
                hl_corrected_list,
                h_tp_list,
                molFracH2O, molFracCO2,
                P_H2O_list, P_CO2_list,
                q_rad_list_CO2, q_rad_list_H2O,
                q_rad_list, q_tot_list,
                CHF_Meyer_list, CHF_Tong_list,

                # Coolant properties
                coolant_velocity_list,
                coolant_reynolds_list,
                coolant_density_list,
                coolant_cond_list,
                coolant_cp_list,
                coolant_visc_list,
                coolant_temp_list,
                coolant_pressure_list,
                coolant_Tsat_list,

                # Wall properties
                hotwall_temp_list,
                coldwall_temp_list,
                wall_cond_list,
                self.wall_material,
                self.wall_thickness,
                hoop_stress_list,
                thermal_stress_list,
                max_wall_stress_list)

            plotter(parameters_plotter, data_plotter)

        #  Writing the results of the study in a CSV file
        if write_in_csv:

            # Write the dimensions of the channels in a CSV file
            file_name = "output/channel_data.csv"
            rows = []
            header = ["Engine z axis", "Engine radius", "Channel width", "Channel height",
                      "Fin width", "Hydraulic diameter", "Channel cross-sectionnal area",
                      "xA", "yA", "xB", "yB", "xC", "yC", "xD", "yD"]
            for i in range(self.nb_points):
                rows.append([
                    z_coord_list[i], r_coord_list[i],
                    channel_width_list[i], channel_height_list[i], effective_fin_thickness_list[i],
                    hydraulic_diameter[i], effective_channel_cross_section[i],
                    channel_vertices['A'][i, 0], channel_vertices['A'][i, 1],
                    channel_vertices['B'][i, 0], channel_vertices['B'][i, 1],
                    channel_vertices['C'][i, 0], channel_vertices['C'][i, 1],
                    channel_vertices['D'][i, 0], channel_vertices['D'][i, 1]
                ])

            with open(file_name, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(header)
                writer.writerows(rows)

            valuexport = open("output/valuexport.csv", "w", newline="")
            valuexport_writer = csv.writer(valuexport)

            # Write header row
            valuexport_writer.writerow([
                "z_coord_list", "r_coord_list", "cross_section_area_list", "effective_fin_thickness_list", "channel_ar_list",
                "channel_width_list", "channel_height_list", "hydraulic_diameter", "effective_channel_cross_section",
                "nb_channels", "gamma_list", "mach_list", "static_pressure_list", "hotgas_total_temp_list",
                "hotgas_recovery_temp_list", "hotgas_static_temp_list", "hotgas_visc_list", "hotgas_cp_list",
                "hotgas_cond_list", "hotgas_pr_list", "hg_list", "sigma_list", "hl_normal_list", "hl_corrected_list",
                "h_tp_list",
                "molFracH2O", "molFracCO2", "P_H2O_list", "P_CO2_list", "q_rad_list_CO2", "q_rad_list_H2O",
                "q_rad_list", "q_tot_list", "coolant_velocity_list", "coolant_reynolds_list", "coolant_density_list",
                "coolant_cond_list", "coolant_cp_list", "coolant_visc_list", "coolant_temp_list", "coolant_pressure_list",
                "coolant_Tsat_list",
                "hotwall_temp_list", "coldwall_temp_list", "wall_cond_list", "CHF_Meyer_list", "CHF_Tong_list", "hoop_stress_list",
                "thermal_stress_list", "max_wall_stress_list"
            ])

            # Collect all rows first
            rows = [
                [
                    z_coord_list[i], r_coord_list[i], cross_section_area_list[i], effective_fin_thickness_list[i], channel_ar_list[i],
                    channel_width_list[i], channel_height_list[i], hydraulic_diameter[i], effective_channel_cross_section[i],
                    self.nb_channels, gamma_list[i], mach_list[i], static_pressure_list[i], hotgas_total_temp_list[i],
                    hotgas_recovery_temp_list[i], hotgas_static_temp_list[i], hotgas_visc_list[i], hotgas_cp_list[i],
                    hotgas_cond_list[i], hotgas_pr_list[i], hg_list[i], sigma_list[i], hl_normal_list[i], hl_corrected_list[i],
                    h_tp_list[i],
                    molFracH2O[i], molFracCO2[i], P_H2O_list[i], P_CO2_list[i], q_rad_list_CO2[i], q_rad_list_H2O[i],
                    q_rad_list[i], q_tot_list[i], coolant_velocity_list[i], coolant_reynolds_list[i], coolant_density_list[i],
                    coolant_cond_list[i], coolant_cp_list[i], coolant_visc_list[i], coolant_temp_list[i], coolant_pressure_list[i],
                    coolant_Tsat_list[i],
                    hotwall_temp_list[i], coldwall_temp_list[i], wall_cond_list[i], CHF_Meyer_list[i], CHF_Tong_list[i],
                    hoop_stress_list[i], thermal_stress_list[i], max_wall_stress_list[i]
                ]
                for i in range(self.nb_points)
            ]
            valuexport_writer.writerows(rows)
            valuexport.close()

        return max_wall_temp, avg_wall_temp, max_stress, coolant_pressure_drop, coolant_temp_increase


tm = TaskManager("input/input_bulk.json")

"""
Created on Sat Jul 12 2025

Original author: Paul M
"""

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
        self.configurations_df = pd.DataFrame(combinations, columns=(self.settings.keys()))
        self.configurations_df["cte_instance"] = pd.Series(dtype="object")
        self.configurations_df["max_wall_temp"] = pd.Series(dtype="float")
        self.configurations_df["avg_wall_temp"] = pd.Series(dtype="float")
        self.configurations_df["max_stress"] = pd.Series(dtype="float")
        self.configurations_df["coolant_pressure_drop"] = pd.Series(dtype="float")
        self.configurations_df["coolant_temp_increase"] = pd.Series(dtype="float")

        # Compute the chamber flow first (it is assumed it is the same for every configuration)
        cte = CoolTheEngine(dict(self.configurations_df.loc[0]))
        chamber_flow_data = cte.compute_chamber_flow()

        # Compute every configuration
        for i, row in self.configurations_df.iterrows():
            # Create and store the instance of CoolTheEngine in the dataframe
            cte = CoolTheEngine(dict(row))
            self.configurations_df.loc[i, "cte_instance"] = cte

            # Compute and unpack results
            try:
                (max_wall_temp, avg_wall_temp, max_stress,
                    coolant_pressure_drop, coolant_temp_increase) = cte.compute_all(chamber_flow_data)
            except ValueError as err:
                print(f"Failed computation for config {i}.\t Error: {err}")
                (max_wall_temp, avg_wall_temp, max_stress,
                    coolant_pressure_drop, coolant_temp_increase) = -1, -1, -1, -1, -1

            # Assign results to respective columns
            self.configurations_df.loc[i, "max_wall_temp"] = max_wall_temp
            self.configurations_df.loc[i, "avg_wall_temp"] = avg_wall_temp
            self.configurations_df.loc[i, "max_stress"] = max_stress
            self.configurations_df.loc[i, "coolant_pressure_drop"] = coolant_pressure_drop
            self.configurations_df.loc[i, "coolant_temp_increase"] = coolant_temp_increase

        print(self.configurations_df)
        self.configurations_df.to_csv("output/bulk_results.csv")

    def single_run(self):
        print("Single run")
        cte = CoolTheEngine(self.settings)
        chamber_flow_data = cte.compute_chamber_flow()
        res = cte.compute_all(chamber_flow_data, show_1D=True, show_2D=True, save_plot=True, write_in_csv=True)
        print(f"Maximum wall temperature : {res[0]:.1f} K")
        print(f"Average wall temperature : {res[1]:.1f} K")
        print(f"Maximum wall stress : {res[2]/1e6:.1f} MPa")
        print(f"Coolant pressure drop : {res[3]/1e5:.2f} bar")
        print(f"Coolant temperature increase : {res[4]:.1f} K")


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

    def compute_chamber_flow(self):
        return t.compute_hotgas_flow(self.contour_file,
                                     self.nb_points,
                                     self.chamber_pressure,
                                     self.ox_mfr,
                                     self.fuel_mfr,
                                     self.ox_name,
                                     self.fuel_name)

    def compute_all(self, chamber_flow_data, show_1D=False, show_2D=False, save_plot=False,
                    write_in_csv=False, mach_list=None):

        # Unpack all the chamber flow data
        (z_coord_list, r_coord_list, throat_diam, throat_area,
            throat_curv_radius, Cstar, Tc, MolWt, hotgas_visc_list, hotgas_cp_list,
            hotgas_cond_list, hotgas_pr_list, cross_section_area_list, gamma_list,
            mach_list, static_pressure_list, molFracCO2, molFracH2O, P_CO2_list, P_H2O_list,
            hotgas_static_temp_list, hotgas_total_temp_list, hotgas_recovery_temp_list,
            x_chamber_throat_exit) = chamber_flow_data

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

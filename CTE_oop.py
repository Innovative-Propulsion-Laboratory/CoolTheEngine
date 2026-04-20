"""
Created on Sat Jul 12 2025

Original author: Paul M
"""

import csv
import json
from pathlib import Path
from typing import Callable
import pandas as pd

# Calculations
import numpy as np
import cte_tools as t
from solver import solver
from itertools import product
import fluid_properties as flp

# Data
from channels import generate_channels

# Graphics
from plotter import plotter
import matplotlib.pyplot as plt



class TaskManager():
    def __init__(
        self,
        json_path: str,
        contour_path: str,
        output_dir: str = "output",
        progress_callback: Callable[[int, int], None] | None = None,
        log_callback: Callable[[str], None] | None = None,
    ):
        with open(json_path, "r") as f:
            self.settings = json.load(f)

        self.contour_path = contour_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.progress_callback = progress_callback
        self.log_callback = log_callback

    def _log(self, message: str) -> None:
        print(message)
        if self.log_callback is not None:
            self.log_callback(message)

    def _report_progress(self, current: int, total: int) -> None:
        if self.progress_callback is not None:
            self.progress_callback(current, total)

    def run(self):
        # Check if any parameter is a list, then run bulk mode
        if any(isinstance(param, list) for param in list(self.settings.values())):
            figs = self.bulk_run()
            bulk_run = True

        # Otherwise, run the single mode (no parameter is a list)
        else:
            figs = self.single_run()
            bulk_run = False

        return bulk_run, figs
    
    def bulk_run(self):
        # Transform into singletons if not a list, and store all parameters in a list of lists
        value_lists = [v if isinstance(v, list) else [v] for v in (self.settings.values())]

        nb_configs = np.prod([len(lst) for lst in value_lists])
        # if nb_configs > 5000:
        #     raise ValueError(f"Too many configurations ({nb_configs}>5000)")

        self._log(f"Solving {nb_configs} unique configurations.")
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

        interp_cache: dict[str, tuple] = {}

        def get_interpolation_functions(coolant_name: str):
            coolant_key = str(coolant_name)
            if coolant_key in interp_cache:
                return interp_cache[coolant_key]

            self._log(f"Generating fluid data tables for {coolant_key}")
            (T_grid, logP_grid, P_grid, Tsat1D, dens_table, cp_table, cond_table,
             visc_table) = flp.build_tables(Tmin=273, Tmax=600, Pmin_bar=1,
                                            Pmax_bar=60, dT=2, nP=100, fluid=coolant_key)

            interp_cache[coolant_key] = flp.make_interpolators(
                T_grid,
                logP_grid,
                P_grid,
                Tsat1D,
                dens_table,
                cp_table,
                cond_table,
                visc_table,
            )
            return interp_cache[coolant_key]

        self._log("Computing chamber flow")
        # Compute the chamber flow first (it is assumed it is the same for every configuration)
        cte = CoolTheEngine(dict(self.configurations_df.loc[0]), self.contour_path)
        chamber_flow_data = cte.compute_chamber_flow()

        # Prepare lists to collect results
        n = len(self.configurations_df)
        cte_instances = [None] * n
        max_wall_temps = [None] * n
        avg_wall_temps = [None] * n
        max_stresses = [None] * n
        coolant_pressure_drops = [None] * n
        coolant_temp_increases = [None] * n

        # Loop only fills Python lists
        total_configs = int(self.configurations_df.shape[0])
        self._report_progress(0, total_configs)
        for i, row in enumerate(self.configurations_df.itertuples(index=False)):
            params = row._asdict()  # convert namedtuple to dict
            cte = CoolTheEngine(params, self.contour_path)
            cte_instances[i] = cte

            try:
                interpolation_functions = get_interpolation_functions(params["coolant_name"])
                res = cte.compute_all(chamber_flow_data, interp_funcs=interpolation_functions)
            except ValueError as err:
                self._log(f"Failed computation for config {i}.\t Error: {err}")
                res = (-1, -1, -1, -1, -1, None)

            (max_wall_temps[i], avg_wall_temps[i], max_stresses[i],
             coolant_pressure_drops[i], coolant_temp_increases[i]) = res
            self._report_progress(i + 1, total_configs)

        self.configurations_df["cte_instance"] = cte_instances
        self.configurations_df["max_wall_temp"] = max_wall_temps
        self.configurations_df["avg_wall_temp"] = avg_wall_temps
        self.configurations_df["max_stress"] = max_stresses
        self.configurations_df["coolant_pressure_drop"] = coolant_pressure_drops
        self.configurations_df["coolant_temp_increase"] = coolant_temp_increases

        fig, axs = plt.subplots(2, 2, dpi=150)
        scatter_data = [
            (coolant_pressure_drops, max_stresses, "Coolant pressure drop [Pa]", "Maximum wall stress [MPa]"),
            (coolant_pressure_drops, avg_wall_temps, "Coolant pressure drop [Pa]", "Average wall temperature [K]"),
            (coolant_pressure_drops, max_wall_temps, "Coolant pressure drop [Pa]", "Maximum wall temperature [K]"),
            (coolant_pressure_drops, coolant_temp_increases, "Coolant pressure drop [Pa]", "Coolant temperature increase [K]")
        ]

        scatters = []
        annots = []
        axes = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]
        for i, (ax, (x, y, xlabel, ylabel)) in enumerate(zip(axes, scatter_data)):
            sc = ax.scatter(x, y, s=1, color="k")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            # Annotation only shows config number
            annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                                bbox=dict(boxstyle="round", fc="w"),
                                arrowprops=dict(arrowstyle="->"), fontsize=8)
            annot.set_visible(False)
            scatters.append(sc)
            annots.append(annot)

        # Add a text box in the top left corner of axs[0,0] for detailed info
        info_text = axs[0, 0].text(0.01, 0.99, "", transform=axs[0, 0].transAxes, va='top', ha='left', fontsize=9,
                                   bbox=dict(boxstyle="round", fc="w", alpha=0.7))
        info_text.set_visible(False)

        def update_all_annots(idx):
            # Update all annotations to the same index
            for sc, annot in zip(scatters, annots):
                pos = sc.get_offsets()[idx]
                annot.xy = pos
                # Only show config number in annotation
                annot.set_text(f"Config {idx}")
                annot.get_bbox_patch().set_facecolor('0.8')
                annot.get_bbox_patch().set_alpha(0.4)

            # Update the info text in the top left of axs[0,0]
            info = (
                f"Config {idx}\n"
                f"Channel width injection : {self.configurations_df.loc[idx, 'channel_widths_inj']*1000:.1f}mm\n"
                f"Channel width converging : {self.configurations_df.loc[idx, 'channel_widths_conv']*1000:.1f}mm\n"
                f"Channel width throat : {self.configurations_df.loc[idx, 'channel_widths_throat']*1000:.1f}mm\n"
                f"Channel width exit : {self.configurations_df.loc[idx, 'channel_widths_exit']*1000:.1f}mm\n"
                f"Channel height injection : {self.configurations_df.loc[idx, 'channel_heights_inj']*1000:.1f}mm\n"
                f"Channel height converging : {self.configurations_df.loc[idx, 'channel_heights_conv']*1000:.1f}mm\n"
                f"Channel height throat : {self.configurations_df.loc[idx, 'channel_heights_throat']*1000:.1f}mm\n"
                f"Channel height exit : {self.configurations_df.loc[idx, 'channel_heights_exit']*1000:.1f}mm\n"
                f"Channel angle injection : {self.configurations_df.loc[idx, 'channel_angles_inj']:.1f}°\n"
                f"Channel angle converging : {self.configurations_df.loc[idx, 'channel_angles_conv']:.1f}°"
            )
            info_text.set_text(info)
            info_text.set_visible(True)

        def set_all_annots_visible(visible):
            for annot in annots:
                annot.set_visible(visible)
            # Also hide/show the info text box
            info_text.set_visible(visible)

        def hover(event):
            found = False
            idx = None
            # Check all axes for a hovered point
            for sc, ax in zip(scatters, axes):
                if event.inaxes == ax:
                    cont, ind = sc.contains(event)
                    if cont:
                        idx = ind["ind"][0]
                        found = True
                        break
            if found and idx is not None:
                update_all_annots(idx)
                set_all_annots_visible(True)
                fig.canvas.draw_idle()
            else:
                set_all_annots_visible(False)
                fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", hover)

        # plt.show()

        self._log(str(self.configurations_df))
        self.configurations_df.to_csv(self.output_dir / "bulk_results.csv", index=False)
    
        return fig

    def single_run(self):
        self._log("Single run")
        cte = CoolTheEngine(self.settings, self.contour_path)
        chamber_flow_data = cte.compute_chamber_flow()
        max_wall_temp, avg_wall_temp, max_stress, coolant_pressure_drop, coolant_temp_increase, figs = cte.compute_all(
            chamber_flow_data,
            show_1D=True,
            show_2D=True,
            save_plot=True,
            write_in_csv=True,
            return_figs=True,
            output_dir=self.output_dir,
        )
        self._log(f"Maximum wall temperature : {max_wall_temp:.1f} K")
        self._log(f"Average wall temperature : {avg_wall_temp:.1f} K")
        self._log(f"Maximum wall stress : {max_stress/1e6:.1f} MPa")
        self._log(f"Coolant pressure drop : {coolant_pressure_drop/1e5:.2f} bar")
        self._log(f"Coolant temperature increase : {coolant_temp_increase:.1f} K")
        
        return figs

class CoolTheEngine():
    def __init__(self, params: dict, contour_path: str):
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
        self.TBC_thickness = float(params["tbc_thickness"])

        # Channel parameters
        self.wall_material = params["wall_material"]
        self.channel_roughness = float(params["channel_roughness"])
        self.nb_channels = int(params["nb_channels"])
        self.wall_thickness = float(params["wall_thickness"])
        self.total_wall_thickness = float(params["total_wall_thickness"])

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
        self.contour_file = contour_path  # Engine contour

    def compute_chamber_flow(self):
        return t.compute_hotgas_flow(self.contour_file,
                                     self.nb_points,
                                     self.chamber_pressure,
                                     self.ox_mfr,
                                     self.fuel_mfr,
                                     self.ox_name,
                                     self.fuel_name)

    def compute_all(self, chamber_flow_data, show_1D=False, show_2D=False, save_plot=False,
                    write_in_csv=False, mach_list=None, interp_funcs=None, return_figs=False,
                    output_dir: str | Path = "output"):

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

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
                        self.coolant_name, self.coolant_mfr, self.TBC_thickness)
        data_channel = (self.nb_channels, channel_width_list,
                        channel_height_list, effective_fin_thickness_list,
                        self.wall_thickness, hydraulic_diameter,
                        effective_channel_cross_section,
                        channel_centerline, beta_list)
        data_chamber = (z_coord_list, r_coord_list, throat_diam, throat_curv_radius, throat_area,
                        self.channel_roughness, cross_section_area_list, mach_list, self.wall_material)

        # Call the main solving loop
        hl_corrected_list, h_tp_list, hg_list, tbc_temp_list, \
            hotwall_temp_list, coldwall_temp_list, q_tot_list, sigma_list, \
            coolant_reynolds_list, coolant_temp_list, coolant_visc_list, \
            coolant_cond_list, coolant_cp_list, coolant_density_list, \
            coolant_velocity_list, coolant_pressure_list, coolant_Tsat_list, wall_cond_list, \
            hl_normal_list, hl_corrected_list, q_rad_list, q_rad_list_CO2, q_rad_list_H2O, \
            CHF_Meyer_list, CHF_Tong_list = solver(data_hotgas, data_coolant, data_channel, data_chamber, interp_funcs)

        hoop_stress_list, thermal_stress_list, max_wall_stress_list = t.compute_1D_wall_stress(self.wall_material, self.wall_thickness, z_coord_list,
                                                                                               r_coord_list, static_pressure_list,
                                                                                               coolant_pressure_list, hotwall_temp_list,
                                                                                               coldwall_temp_list, self.nb_channels,
                                                                                               effective_fin_thickness_list, channel_height_list,
                                                                                               channel_width_list, self.total_wall_thickness)

        max_wall_temp = np.max(hotwall_temp_list)
        avg_wall_temp = np.average(hotwall_temp_list)
        max_stress = np.max(max_wall_stress_list)
        coolant_pressure_drop = coolant_pressure_list[-1] - coolant_pressure_list[0]
        coolant_temp_increase = coolant_temp_list[0] - coolant_temp_list[-1]

        # PLot the results
        if show_1D or show_2D:
            # Display of the 1D analysis results
            parameters_plotter = (show_1D, show_2D, self.figure_dpi, save_plot, self)

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
                tbc_temp_list,
                hotwall_temp_list,
                coldwall_temp_list,
                wall_cond_list,
                self.wall_material,
                self.wall_thickness,
                hoop_stress_list,
                thermal_stress_list,
                max_wall_stress_list)

            figs = plotter(parameters_plotter, data_plotter, output_dir=output_dir)

        #  Writing the results of the study in a CSV file
        if write_in_csv:

            # Write the dimensions of the channels in a CSV file
            file_name = output_dir / "channel_data.csv"
            rows = []
            header = ["Engine z axis", "Engine radius", "Channel width", "Channel height",
                      "Fin width", "Hydraulic diameter", "Channel cross-sectionnal area",
                      "xA", "yA", "xB", "yB", "xC", "yC", "xD", "yD"]
            for i in range(self.nb_points):
                rows.append([
                    z_coord_list[i]*1000, r_coord_list[i]*1000,
                    channel_width_list[i]*1000, channel_height_list[i]*1000, effective_fin_thickness_list[i]*1000,
                    hydraulic_diameter[i]*1000, effective_channel_cross_section[i],
                    channel_vertices['A'][i, 0]*1000, channel_vertices['A'][i, 1]*1000,
                    channel_vertices['B'][i, 0]*1000, channel_vertices['B'][i, 1]*1000,
                    channel_vertices['C'][i, 0]*1000, channel_vertices['C'][i, 1]*1000,
                    channel_vertices['D'][i, 0]*1000, channel_vertices['D'][i, 1]*1000
                ])

            with open(file_name, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(header)
                writer.writerows(rows)

            valuexport = open(output_dir / "valuexport.csv", "w", newline="")
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

        if return_figs:
            return max_wall_temp, avg_wall_temp, max_stress, coolant_pressure_drop, coolant_temp_increase, figs
        else:
            return max_wall_temp, avg_wall_temp, max_stress, coolant_pressure_drop, coolant_temp_increase


if __name__ == "__main__":
    tm = TaskManager("input/input.json", "input/redwing_contour.csv", output_dir="output")
    bulk_run, figs = tm.run()
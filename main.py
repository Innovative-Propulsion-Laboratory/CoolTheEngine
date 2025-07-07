"""
Created on Fri Nov 27 14:47:27 2020

Original author: Julien S (2020)

Refactored and improved by Mehdi D, Paul B, Paul M, Eve X and Antoine R (2023)

Improved and enhanced by Paul Marchi (july 2025)
"""

from scipy.interpolate import interp1d
import time
import csv
import json

# Calculations
import numpy as np
import cte_tools as t
from solver import solver

# Data
from channels import generate_channels
import cea

# Graphics
import matplotlib.pyplot as plt
from tqdm import tqdm  # For progress bars
from plotter import plotter


start_time = time.perf_counter()  # Beginning of the timer

print("██████████████████████████  Cool The Engine  V3.0  █████████████████████████")
print("█                                                                          █")
print("█                  Innovative Propulsion Laboratory - IPL                  █")
print("█__________________________________________________________________________█")
print("█                                                                          █")
print("█ Initialisation                                                           █")
print("█                                                                          █")

# Initial definitions

input_file = "input/input.txt"  # Engine parameters
contour_file = "input/engine_contour.csv"  # Engine contour

# Constant input_data_list
figure_dpi = 150  # Dots Per Inch (DPI) for all figures (lower=faster)
plot_detail = 3  # 0=No plots; 1=Important plots; 2=Less important plots: 3=All plots
show_3d_plots = False
show_2D_temperature = True
do_final_3d_plot = False
write_in_csv = True
nb_points = 100
lim_3D_plot = 0.05
save_plots = True

# %% Reading input data
contour_data = np.genfromtxt(contour_file, delimiter=",", skip_header=1)
z_coord_list = contour_data[0:, 0]/1000
r_coord_list = contour_data[0:, 1]/1000
nb_points_raw = len(z_coord_list)  # Number of points
nb_points = 500

# Reduce the number of points using scipy 1D interpolation

# Create new x values evenly spaced between the min and max of the original x_coord_list
x_new = np.linspace(z_coord_list[0], z_coord_list[-1], nb_points)

# Interpolate y values at the new x positions
interp_func = interp1d(z_coord_list, r_coord_list, kind='linear')
y_new = interp_func(x_new)

# Replace the original lists with the interpolated ones
z_coord_list = x_new
r_coord_list = y_new

# Read input data from JSON file
with open("input/input.json", "r") as f:
    data = json.load(f)

# Engine parameters
chamber_pressure = float(data["chamber_pressure"])
coolant_inlet_pressure = float(data["P_coolant"])
coolant_inlet_temp = float(data["T_coolant"])
ox_mfr = float(data["ox_mfr"])
fuel_mfr = float(data["fuel_mfr"])
coolant_mfr = float(data["coolant_mfr"])
ox_name = data["oxidizer_name"]
fuel_name = data["fuel_name"]
coolant_name = data["coolant_name"]

# Channel parameters
wall_material = data["wall_material"]
channel_roughness = float(data["channel_roughness"])
nb_channels = int(data["nb_channels"])
wall_thickness = float(data["wall_thickness"])

# Widths
width_inj = float(data["channel_widths"]["inj"])
width_conv = float(data["channel_widths"]["conv"])
width_throat = float(data["channel_widths"]["throat"])
width_exit = float(data["channel_widths"]["exit"])

# Heights
ht_inj = float(data["channel_heights"]["inj"])
ht_conv = float(data["channel_heights"]["conv"])
ht_throat = float(data["channel_heights"]["throat"])
ht_exit = float(data["channel_heights"]["exit"])

# Angles
beta_inj = float(data["channel_angles"]["inj"])
beta_conv = float(data["channel_angles"]["conv"])
beta_throat = float(data["channel_angles"]["throat"])
beta_exit = float(data["channel_angles"]["exit"])

chamber_radius = r_coord_list[0]  # Radius of the chamber (in m)
throat_radius = np.min(r_coord_list)  # Radius of the throat (in m)
exit_radius = r_coord_list[-1]  # Radius of the exit (in m)
throat_diam = 2 * throat_radius  # Diameter of the throat (in m)
chamber_diam = 2 * chamber_radius  # Diameter of the chamber (in m)
exit_diam = 2 * exit_radius  # Diameter of the exit (in m)

chamber_area = np.pi * chamber_radius**2  # Area of the chamber (in m²)
throat_area = np.pi * throat_radius**2  # Area of the throat (in m²)
exit_area = np.pi * exit_radius**2  # Area of the exit (in m²)

contraction_ratio = chamber_area / throat_area  # Contraction ratio (A/A*)
expansion_ratio = exit_area / throat_area  # Expansion ratio (Ae/A*)
throat_curv_radius = 1.5 * throat_radius  # Curvature radius before the throat (in m)

i_throat = np.argmin(np.abs(r_coord_list))
i_convergent = 0
for i in range(1, len(r_coord_list)):
    if r_coord_list[i] < r_coord_list[i-1]:
        i_convergent = i
        break

Cstar, Tc, MolWt = cea.compute_Cstar_Tc_MolWt(chamber_pressure, ox_mfr/fuel_mfr,
                                              ox_name, fuel_name, expansion_ratio)

hotgas_mu_chamber, hotgas_cp_chamber, hotgas_lambda_chamber, hotgas_pr_chamber, \
    hotgas_mu_throat, hotgas_cp_throat, hotgas_lambda_throat, hotgas_pr_throat, \
    hotgas_mu_exit, hotgas_cp_exit, hotgas_lambda_exit, hotgas_pr_exit\
    = cea.get_hotgas_properties(chamber_pressure, ox_mfr/fuel_mfr,
                                ox_name, fuel_name, expansion_ratio)

x_chamber_throat_exit = [z_coord_list[0], z_coord_list[i_convergent],
                         z_coord_list[i_throat], z_coord_list[-1]]

hotgas_visc_list = interp1d(x_chamber_throat_exit,
                            [hotgas_mu_chamber, hotgas_mu_chamber, hotgas_mu_throat, hotgas_mu_exit],
                            kind='linear')(z_coord_list)
hotgas_cp_list = interp1d(x_chamber_throat_exit,
                          [hotgas_cp_chamber, hotgas_cp_chamber, hotgas_cp_throat, hotgas_cp_exit],
                          kind='linear')(z_coord_list)
hotgas_cond_list = interp1d(x_chamber_throat_exit,
                            [hotgas_lambda_chamber, hotgas_lambda_chamber, hotgas_lambda_throat, hotgas_lambda_exit],
                            kind='linear')(z_coord_list)
hotgas_pr_list = interp1d(x_chamber_throat_exit,
                          [hotgas_pr_chamber, hotgas_pr_chamber, hotgas_pr_throat, hotgas_pr_exit],
                          kind='linear')(z_coord_list)


# Computation of the cross-sectional area along the engine
cross_section_area_list = np.pi*r_coord_list**2

# Adiabatic constant (gamma) parametrization
print("█ Computing gamma                                                          █")

gamma_list = cea.compute_gamma(chamber_pressure, ox_mfr/fuel_mfr,
                               ox_name, fuel_name, cross_section_area_list/throat_area)

# Computation of mach number of the hot gases
mach_list = cea.compute_mach(chamber_pressure, ox_mfr/fuel_mfr,
                             ox_name, fuel_name, cross_section_area_list/throat_area)

# Static pressure computation
static_pressure_list = np.zeros_like(z_coord_list)  # (in Pa)

with tqdm(total=nb_points,
          desc="█ Computing static pressure    ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points):
        static_pressure_list[i] = t.pressure_solv(mach_list[i], gamma_list[i], chamber_pressure)
        progressbar.update(1)

# Partial pressure computation and interpolation of the molar fraction
molFracH2O_chamber, molFracH2O_throat, molFracH2O_exit\
    = cea.compute_H2O_molar_fractions(chamber_pressure, ox_mfr/fuel_mfr,
                                      ox_name, fuel_name, expansion_ratio)
molFracCO2_chamber, molFracCO2_throat, molFracCO2_exit\
    = cea.compute_CO2_molar_fractions(chamber_pressure, ox_mfr/fuel_mfr,
                                      ox_name, fuel_name, expansion_ratio)

# Linear interpolation of molar fractions for H2O and CO2
molFracH2O = interp1d(x_chamber_throat_exit, [molFracH2O_chamber, molFracH2O_chamber, molFracH2O_throat, molFracH2O_exit], kind='linear')(z_coord_list)
molFracCO2 = interp1d(x_chamber_throat_exit, [molFracCO2_chamber, molFracCO2_chamber, molFracCO2_throat, molFracCO2_exit], kind='linear')(z_coord_list)

# Partial pressure of H2O and CO2
P_H2O_list = np.array([static_pressure_list[i] * molFracH2O[i] for i in range(0, nb_points)])
P_CO2_list = np.array([static_pressure_list[i] * molFracCO2[i] for i in range(0, nb_points)])

# Hot gas temperature computation
hotgas_static_temp_list = np.zeros_like(z_coord_list)  # List of static hot gas temperatures
with tqdm(total=nb_points,
          desc="█ Computing gas static_temperature    ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points):
        hotgas_static_temp_list[i] = t.temperature_hotgas_solv(mach_list[i], gamma_list[i], Tc)
        progressbar.update(1)

    # We assume that the total temperature is constant
    hotgas_total_temp_list = Tc*np.ones_like(z_coord_list)  # List of total hot gas temperatures

# Computation of adiabatic wall temperature (recovery temperature)
hotgas_recovery_temp_list = np.array([t.get_recovery_temperature(hotgas_total_temp_list[i], gamma_list[i], mach_list[i], hotgas_pr_list[i]) for i in
                                      range(0, nb_points)])

# Dimensions
print("█ Computing channel geometry                                               █")
print("█                                                                          █")

# Pack the data in tuples
profile = (z_coord_list, r_coord_list)
widths = (width_inj, width_conv, width_throat, width_exit)
heights = (ht_inj, ht_conv, ht_throat, ht_exit)
angles = (beta_inj, beta_conv, beta_throat, beta_exit)

# Generate the cooling channels
channel_vertices, channel_centerline, channel_inclination, channel_width_list, channel_height_list, \
    initial_channel_cross_section, effective_channel_cross_section, \
    hydraulic_diameter, initial_fin_thickness, effective_fin_thickness_list, \
    alpha_list, beta_list, channel_total_length = generate_channels(profile, widths, heights,
                                                                    angles, wall_thickness,
                                                                    nb_channels, x_chamber_throat_exit)

# Write the dimensions of the channels in a CSV file
file_name = "output/channel_data.csv"
with open(file_name, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(("Engine z axis", "Engine radius", "Channel width", "Channel height",
                     "Fin width", "Hydraulic diameter", "Channel cross-sectionnal area,"
                     "xA", "yA", "xB", "yB", "xC", "yC", "xD", "yD",))
    for i in range(nb_points):
        writer.writerow((z_coord_list[i], r_coord_list[i],
                         channel_width_list[i], channel_height_list[i], effective_fin_thickness_list[i],
                         hydraulic_diameter[i], effective_channel_cross_section[i],
                         channel_vertices['A'][:, 0], channel_vertices['A'][:, 1],
                         channel_vertices['B'][:, 0], channel_vertices['B'][:, 1],
                         channel_vertices['C'][:, 0], channel_vertices['C'][:, 1],
                         channel_vertices['D'][:, 0], channel_vertices['D'][:, 1]))

end_init_time = time.perf_counter()  # End of the initialisation timer
time_elapsed = f"{round(end_init_time - start_time, 2)}"  # Initialisation elapsed time (in s)
time_elapsed_i = f"{time_elapsed} s"

start_main_time = time.perf_counter()  # Start of the main solution timer


# Main computation

data_hotgas = (hotgas_recovery_temp_list, hotgas_static_temp_list, hotgas_visc_list, hotgas_pr_list,
               hotgas_cp_list, hotgas_cond_list, MolWt, gamma_list, chamber_pressure,
               Cstar, P_H2O_list, P_CO2_list)
data_coolant = (coolant_inlet_temp, coolant_inlet_pressure, coolant_name, coolant_mfr)
data_channel = (nb_channels, channel_width_list, channel_height_list, effective_fin_thickness_list,
                wall_thickness, hydraulic_diameter, effective_channel_cross_section,
                channel_centerline, beta_list)
data_chamber = (z_coord_list, r_coord_list, throat_diam, throat_curv_radius, throat_area,
                channel_roughness, cross_section_area_list, mach_list, wall_material)

# Call the main solving loop
hl_corrected_list, hg_list, \
    hotwall_temp_list, coldwall_temp_list, q_tot_list, sigma_list, \
    coolant_reynolds_list, coolant_temp_list, coolant_visc_list, \
    coolant_cond_list, coolant_cp_list, coolant_density_list, \
    coolant_velocity_list, coolant_pressure_list, wall_cond_list, hg_list, \
    hl_normal_list, hl_corrected_list, q_rad_list, q_rad_list_CO2, q_rad_list_H2O \
    = solver(data_hotgas, data_coolant, data_channel, data_chamber)

end_m = time.perf_counter()  # End of the main solution timer
time_elapsed = f"{round(end_m - start_main_time, 2)}"  # Main computation elapsed time (in s)
if len(time_elapsed) <= 3:
    time_elapsed_m = f"   {time_elapsed} s"
elif len(time_elapsed) == 4:
    time_elapsed_m = f"  {time_elapsed} s"
elif len(time_elapsed) == 5:
    time_elapsed_m = f" {time_elapsed} s"
else:
    time_elapsed_m = f"{time_elapsed} s"

# Display of the 1D analysis results
print("█                                                                          █")


start_p = time.perf_counter()  # Start of the end timer
parameters_plotter = (plot_detail, show_2D_temperature, do_final_3d_plot,
                      figure_dpi, True, save_plots)


data_plotter = (  # Engine geometry
    z_coord_list * 1000,
    r_coord_list * 1000,
    z_coord_list,
    r_coord_list,
    cross_section_area_list,

    # Channel geometry
    initial_fin_thickness,
    effective_fin_thickness_list,
    channel_width_list,
    channel_height_list,
    hydraulic_diameter,
    initial_channel_cross_section,
    effective_channel_cross_section,
    nb_channels,
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
    molFracH2O, molFracCO2,
    P_H2O_list, P_CO2_list,
    q_rad_list_CO2, q_rad_list_H2O,
    q_rad_list, q_tot_list,

    # Coolant properties
    coolant_velocity_list,
    coolant_reynolds_list,
    coolant_density_list,
    coolant_cond_list,
    coolant_cp_list,
    coolant_visc_list,
    coolant_temp_list,
    coolant_pressure_list,

    # Wall properties
    hotwall_temp_list,
    coldwall_temp_list,
    wall_cond_list,
    wall_material,
    wall_thickness)

plotter(parameters_plotter, data_plotter)
end_p = time.perf_counter()
time_elapsed = f"{round(end_p - start_p, 2)}"  # Total elapsed time
if len(time_elapsed) <= 3:
    time_elapsed_p = f"   {time_elapsed} s"
elif len(time_elapsed) == 4:
    time_elapsed_p = f"  {time_elapsed} s"
elif len(time_elapsed) == 5:
    time_elapsed_p = f" {time_elapsed} s"
else:
    time_elapsed_p = f"{time_elapsed} s"


start_e = time.perf_counter()  # Start of the end timer
#  Writing the results of the study in a CSV file
if write_in_csv:
    print("█ Writing results in .csv files                                            █")
    valuexport = open("output/valuexport.csv", "w", newline="")
    valuexport_writer = csv.writer(valuexport)

    # Write header row
    valuexport_writer.writerow([
        "z_coord_list", "r_coord_list", "cross_section_area_list", "effective_fin_thickness_list",
        "channel_width_list", "channel_height_list", "hydraulic_diameter", "effective_channel_cross_section",
        "nb_channels", "gamma_list", "mach_list", "static_pressure_list", "hotgas_total_temp_list",
        "hotgas_recovery_temp_list", "hotgas_static_temp_list", "hotgas_visc_list", "hotgas_cp_list",
        "hotgas_cond_list", "hotgas_pr_list", "hg_list", "sigma_list", "hl_normal_list", "hl_corrected_list",
        "molFracH2O", "molFracCO2", "P_H2O_list", "P_CO2_list", "q_rad_list_CO2", "q_rad_list_H2O",
        "q_rad_list", "q_tot_list", "coolant_velocity_list", "coolant_reynolds_list", "coolant_density_list",
        "coolant_cond_list", "coolant_cp_list", "coolant_visc_list", "coolant_temp_list", "coolant_pressure_list",
        "hotwall_temp_list", "coldwall_temp_list", "wall_cond_list", "wall_material"
    ])

    for i in range(0, nb_points):
        valuexport_writer.writerow([
            z_coord_list[i], r_coord_list[i], cross_section_area_list[i], effective_fin_thickness_list[i],
            channel_width_list[i], channel_height_list[i], hydraulic_diameter[i], effective_channel_cross_section[i],
            nb_channels, gamma_list[i], mach_list[i], static_pressure_list[i], hotgas_total_temp_list[i],
            hotgas_recovery_temp_list[i], hotgas_static_temp_list[i], hotgas_visc_list[i], hotgas_cp_list[i],
            hotgas_cond_list[i], hotgas_pr_list[i], hg_list[i], sigma_list[i], hl_normal_list[i], hl_corrected_list[i],
            molFracH2O[i], molFracCO2[i], P_H2O_list[i], P_CO2_list[i], q_rad_list_CO2[i], q_rad_list_H2O[i],
            q_rad_list[i], q_tot_list[i], coolant_velocity_list[i], coolant_reynolds_list[i], coolant_density_list[i],
            coolant_cond_list[i], coolant_cp_list[i], coolant_visc_list[i], coolant_temp_list[i], coolant_pressure_list[i],
            hotwall_temp_list[i], coldwall_temp_list[i], wall_cond_list[i], wall_material
        ])

    valuexport.close()

# Execution time display
end_t = time.perf_counter()  # End of the total timer
time_elapsed = f"{round(end_t - start_e, 2)}"  # End elapsed time (in s)
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

print("█                                                                          █")
print("█__________________________________________________________________________█")
print("█                                                                          █")
print(f"█ Execution time for the initialisation       : {time_elapsed_i}                     █")
print("█                                                                          █")
print(f"█ Execution time for the main computation     : {time_elapsed_m}                   █")

if plot_detail >= 1:
    print("█                                                                          █")
    print(f"█ Execution time for the plotting             : {time_elapsed_p}                   █")

print("█                                                                          █")
print(f"█ Execution time for the end of the program   : {time_elapsed_e}                   █")
print("█                                                                          █")
print(f"█ Total execution time                        : {time_elapsed_t}                   █")
print("█                                                                          █")
print("███████████████████████████████████ END ████████████████████████████████████")

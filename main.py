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
from solver2D import carto2D
from volume3d import carto3d, view3d
import matplotlib.pyplot as plt
from tqdm import tqdm  # For progress bars

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

# %% Reading input data
contour_data = np.genfromtxt(contour_file, delimiter=",", skip_header=1)
z_coord_list = contour_data[0:, 0]/1000
r_coord_list = contour_data[0:, 1]/1000
nb_points_raw = len(z_coord_list)  # Number of points
nb_points = 100

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
throat_curv_radius = 1.5 * chamber_radius  # Curvature radius before the throat (in m)

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
hotgas_conv_list = interp1d(x_chamber_throat_exit,
                            [hotgas_lambda_chamber, hotgas_lambda_chamber, hotgas_lambda_throat, hotgas_lambda_exit],
                            kind='linear')(z_coord_list)
hotgas_pr_list = interp1d(x_chamber_throat_exit,
                          [hotgas_pr_chamber, hotgas_pr_chamber, hotgas_pr_throat, hotgas_pr_exit],
                          kind='linear')(z_coord_list)

# Plot of the engine profile
if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, r_coord_list, color='black')
    plt.title('Profile of the engine (left : chamber and right : divergent)', color='black')
    plt.show()

# Computation of the cross-sectional area along the engine
cross_section_area_list = np.pi*r_coord_list**2

# Plot of the cross-sectional area in the engine
if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, cross_section_area_list, color='black')
    plt.title("Cross section area of engine (in m²) as a function of engine axis")
    plt.show()

# Adiabatic constant (gamma) parametrization
print("█ Computing gamma                                                          █")

gamma_list = cea.compute_gamma(chamber_pressure, ox_mfr/fuel_mfr,
                               ox_name, fuel_name, cross_section_area_list/throat_area)

# Plot of the gamma linearisation
if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, gamma_list, color='gold')
    plt.title("Gamma of hot gases as a function of engine axis")
    plt.show()

# Computation of mach number of the hot gases
mach_list = cea.compute_mach(chamber_pressure, ox_mfr/fuel_mfr,
                             ox_name, fuel_name, cross_section_area_list/throat_area)

# Plots of the Mach number in the engine (2D/3D)
if plot_detail >= 1:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, mach_list, color='gold')
    plt.title("Mach number as a function of the engine axis")
    plt.show()

# Static pressure computation
static_pressure_list = np.zeros_like(z_coord_list)  # (in Pa)

with tqdm(total=nb_points,
          desc="█ Computing static pressure    ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points):
        static_pressure_list[i] = t.pressure_solv(mach_list[i], gamma_list[i], chamber_pressure)
        progressbar.update(1)

# Plot of the static pressure along the engine
if plot_detail >= 2:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, static_pressure_list, color='gold')
    plt.title("Global static pressure (in Pa) as a function of the engine axis")
    plt.show()

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
PH2O_list = np.array([static_pressure_list[i] * molFracH2O[i] for i in range(0, nb_points)])
PCO2_list = np.array([static_pressure_list[i] * molFracCO2[i] for i in range(0, nb_points)])

# Plots of molar fraction and partial pressure
if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, molFracH2O, color='blue', label='H20')
    plt.plot(z_coord_list, molFracCO2, color='orange', label='C02')
    plt.title("Molar fraction of as a function of the engine axis")
    plt.legend(loc='center left')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, PH2O_list, color='blue', label='H20')
    plt.plot(z_coord_list, PCO2_list, color='orange', label='C02')
    plt.title("Partial static pressure (in Pa) of as a function of the engine axis")
    plt.legend(loc='center left')
    plt.show()

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

# Plot of the temperatures in the engine
if plot_detail >= 2:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, hotgas_static_temp_list, label="Static temperature")
    plt.plot(z_coord_list, hotgas_total_temp_list, label="Total temperature")
    plt.plot(z_coord_list, hotgas_recovery_temp_list, label="Static temperature")
    plt.title("Gas temperature (in K) as a function of the engine axis")
    plt.legend()
    plt.show()

# Dimensions
print("█ Computing channel geometry                                               █")
print("█                                                                          █")

# Pack the data in tuples
profile = (z_coord_list, r_coord_list)
widths = (width_inj, width_conv, width_throat, width_exit)
heights = (ht_inj, ht_conv, ht_throat, ht_exit)
angles = (beta_inj, beta_conv, beta_throat, beta_exit)

# Generate the cooling channels
channel_vertices, channel_centerline, channel_inclination, width_list, ht_list, \
    effective_channel_cross_section, hydraulic_diameter, effective_fin_thickness, \
    alpha_list, beta_list, channel_total_length = generate_channels(profile, widths, heights,
                                                                    angles, wall_thickness,
                                                                    nb_channels, x_chamber_throat_exit,
                                                                    plot_detail, write_in_csv,
                                                                    figure_dpi, plot_dir=None)

# Write the dimensions of the channels in a CSV file
file_name = "output/channel_data.csv"
with open(file_name, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(("Engine z axis", "Engine radius", "Channel width", "Channel height",
                     "Fin width", "Hydraulic diameter", "Channel cross-sectionnal area"))
    for i in range(nb_points):
        writer.writerow((z_coord_list[i], r_coord_list[i],
                         width_list[i], ht_list[i], effective_fin_thickness[i],
                         hydraulic_diameter[i], effective_channel_cross_section[i]))

end_init_time = time.perf_counter()  # End of the initialisation timer
time_elapsed = f"{round(end_init_time - start_time, 2)}"  # Initialisation elapsed time (in s)
time_elapsed_i = f"{time_elapsed} s"

start_main_time = time.perf_counter()  # Start of the main solution timer


# Main computation

data_hotgas = (hotgas_recovery_temp_list, hotgas_static_temp_list, hotgas_visc_list, hotgas_pr_list,
               hotgas_cp_list, hotgas_conv_list, MolWt, gamma_list, chamber_pressure,
               Cstar, PH2O_list, PCO2_list)
data_coolant = (coolant_inlet_temp, coolant_inlet_pressure, coolant_name, coolant_mfr)
data_channel = (nb_channels, width_list, ht_list, effective_fin_thickness,
                wall_thickness, hydraulic_diameter, effective_channel_cross_section,
                channel_centerline, beta_list)
data_chamber = (z_coord_list, r_coord_list, throat_diam, throat_curv_radius, throat_area,
                channel_roughness, cross_section_area_list, mach_list, wall_material)

# Call the main solving loop
hl_corrected_list, hg_list, \
    hotwall_temp_list, coldwall_temp_list, q_conv_list, sigma_list, \
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

if plot_detail >= 1:
    start_d1 = time.perf_counter()  # Start of the display of 1D timer
    print("█ Display of results                                                       █")
    print("█                                                                          █")
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, hl_corrected_list, color='blue', label='Hl corrected (Luka Denies)')
    plt.plot(z_coord_list, hl_normal_list, color='cyan', label='Hl')
    plt.title("Convection coeff as a function of the engine axis")
    plt.legend(loc='upper left')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, coldwall_temp_list, color='blue', label='Twl')
    plt.plot(z_coord_list, hotwall_temp_list, color='red', label='Twg')
    plt.title('Wall temperature (in K) as a function of engine axis')
    plt.legend(loc='lower left')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, coolant_temp_list, color='blue')
    plt.title('Coolant temperature (in K) as a function of engine axis')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, q_conv_list, color='red')
    plt.title('Heat flux (in W) as a function of engine axis')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, coolant_pressure_list, color='orange')
    plt.title('Pressure drop in the cooling channels')
    plt.show()

if plot_detail >= 2:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, wall_cond_list, color='orangered')
    plt.title('Conductivity of the wall as a function of engine axis')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, hg_list, color='orangered')
    plt.title('Convection coefficient Hg as a function of engine axis')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, coolant_density_list, color='blue')
    plt.title('Volumic mass of the coolant as a function of engine axis')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, q_rad_list_CO2, color='r', label='CO2')
    plt.plot(z_coord_list, q_rad_list_H2O, color='b', label='H2O')
    plt.plot(z_coord_list, q_rad_list, color='g', label='total')
    plt.title('Radiative heat flux(W/m2)')
    plt.legend()
    plt.show()

if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, coolant_reynolds_list, color='blue')
    plt.title("Reynolds number of the coolant as a function of the engine axis")
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, hotgas_visc_list, color='orangered')
    plt.title("Gas viscosity as a function of the engine axis")
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, hotgas_cp_list, color='orangered')
    plt.title("Gas Cp as a function of the engine axis")
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, hotgas_conv_list, color='orangered')
    plt.title("Gas conductivity as a function of engine axis")
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, hotgas_pr_list, color='orangered')
    plt.title("Gas Prandtl number as a function of engine axis")
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, sigma_list, color='orangered')
    plt.title("Sigma as a function of the engine axis")
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, coolant_cond_list, color='blue')
    plt.title('Conductivity of the coolant as a function of engine axis')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, coolant_cp_list, color='blue')
    plt.title('Cp of the coolant as a function of engine axis')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, coolant_visc_list, color='blue')
    plt.title('Viscosity of the coolant as a function of engine axis')
    plt.show()

if plot_detail >= 1:
    end_d1 = time.perf_counter()  # End of the display of 1D timer
    time_elapsed = f"{round(end_d1 - start_d1, 2)}"  # 1D display elapsed time (in s)
    if len(time_elapsed) <= 3:
        time_elapsed_d1 = f"   {time_elapsed} s"
    elif len(time_elapsed) == 4:
        time_elapsed_d1 = f"  {time_elapsed} s"
    elif len(time_elapsed) == 5:
        time_elapsed_d1 = f" {time_elapsed} s"
    else:
        time_elapsed_d1 = f"{time_elapsed} s"

# Flux computation in 2D and 3D
"""2D flux computation"""

if show_2D_temperature:
    start_d2 = time.perf_counter()  # Start of the display of 2D timer
    # At the beginning of the chamber
    print("█ Results at the beginning of the chamber :                                █")
    dx = 0.00004  # *3.5
    location = " at the beginning of the chamber"
    carto2D(eff[-1] + larg_canal[-1], larg_canal[-1], wall_thickness, ht_canal[-1], dx, hg_list[-1],
            wallcond_list[-1], hotgas_recovery_temp_list[-1], hlcor_list[-1], tempcoolant_list[-1], 5, True, 1,
            location,
            False)

    # At the throat
    print("█ Results at the throat :                                                  █")
    pos_col = ycanaux.index(min(ycanaux))
    dx = 0.000025  # *3.5
    location = " at the throat"
    carto2D(larg_ailette_list[pos_col] + larg_canal[pos_col], larg_canal[pos_col], e_col, ht_canal[pos_col],
            dx, hg_list[pos_col], wallcond_list[pos_col], hotgas_recovery_temp_list[pos_col], hlcor_list[pos_col],
            tempcoolant_list[pos_col], 15, True, 2, location, False)
    # At the end of the divergent
    print("█ Results at the manifold :                                                █")
    dx = 0.00004
    location = " at the manifold"
    carto2D(larg_ailette_list[0] + larg_canal[0], larg_canal[0], e_tore, ht_canal[0], dx, hg_list[0],
            wallcond_list[0], hotgas_recovery_temp_list[0], hlcor_list[0], tempcoolant_list[0], 5, True, 1, location,
            False)

    end_d2 = time.perf_counter()  # End of the display of 2D timer
    time_elapsed = f"{round(end_d2 - start_d2, 2)}"  # 2D display elapsed time (in s)
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
    with tqdm(total=nb_points_channel,
              desc="█ 3D graph computation         ",
              unit="|   █", bar_format="{l_bar}{bar}{unit}",
              ncols=76) as progressbar:
        for i in range(0, nb_points_channel):
            temperature_slice = carto2D(larg_ailette_list[i] + larg_canal[i], larg_canal[i], wall_thickness[i],
                                        ht_canal[i], dx, hg_list[i], wallcond_list[i], hotgas_recovery_temp_list[i],
                                        hlnormal_list[i], tempcoolant_list[i], 3, False, 1, "", True)
            temperature_slice_list.append(temperature_slice)
            progressbar.update(1)

    # Stack all these slices in a final 3D plot
    carto3d([0, 0, 0], z_coord_list, ycanaux, temperature_slice_list, plt.cm.Spectral_r,
            '3D view of wall temperatures (in K)', nb_channels, limitation)
    print("█                                                                          █")
    # End of the 3D display timer
    end_d3 = time.perf_counter()
    time_elapsed = f"{round(end_d3 - start_d3, 2)}"  # 3D display elapsed time (in s)
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
hotgas_recovery_temp_list.reverse()
z_coord_list.reverse()
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
newxhtre = [z_coord_list[0]]
newyhtre = [ycanaux[0] + ht_canal[0]]
for i in range(1, nb_points_channel):
    if i == (nb_points_channel - 1):
        angle = angles[i - 1]
        angles.append(angle)
    else:
        vect1 = (z_coord_list[i] - z_coord_list[i - 1]) / (
                (((ycanaux[i] - ycanaux[i - 1]) ** 2) + ((z_coord_list[i] - z_coord_list[i - 1]) ** 2)) ** 0.5)
        vect2 = (z_coord_list[i + 1] - z_coord_list[i]) / (
                (((ycanaux[i + 1] - ycanaux[i]) ** 2) + ((z_coord_list[i + 1] - z_coord_list[i]) ** 2)) ** 0.5)
        angle1 = np.rad2deg(np.arccos(vect1))
        angle2 = np.rad2deg(np.arccos(vect2))
        angle = angle2
        angles.append(angle)
    newx = z_coord_list[i] + ht_canal[i] * np.sin(np.deg2rad(angles[i]))
    newy = ycanaux[i] + ht_canal[i] * np.cos(np.deg2rad(angles[i]))
    newxhtre.append(newx)
    newyhtre.append(newy)

# Checking the height of channels
verification = []
print("█ Checking and computing channel height                                    █")
for i in range(0, nb_points_channel):
    verifhtre = (((newxhtre[i] - z_coord_list[i]) ** 2) + ((newyhtre[i] - ycanaux[i]) ** 2)) ** 0.5
    verification.append(verifhtre)

if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(newxhtre, newyhtre, color='blue', label='New height')
    plt.plot(z_coord_list, ycanaux, color='chocolate', label='Former height')
    plt.title("Geometrical aspect of the channel (height as a function of the engine axis)")
    plt.axis("equal")
    plt.legend(loc='upper left')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(z_coord_list, verification)
    plt.title("Checking the height of the generated channels")
    plt.show()

# %% Writing the results of the study in a CSV file

if write_in_csv:
    print("█ Writing results in .csv files                                            █")
    valuexport = open("output/valuexport.csv", "w", newline="")
    geometry1 = open("output/geometry1.csv", "w", newline="")
    geometry2 = open("output/geometry2.csv", "w", newline="")
    valuexport_writer = csv.writer(valuexport)
    geometry1_writer = csv.writer(geometry1)
    geometry2_writer = csv.writer(geometry2)
    valuexport_writer.writerow(
        ("Engine x axix", "Engine diameter", "Area of gas engine", "Gas gamma",
         "Mach number", "Gas pressure", "Gas recovery temperature",
         "Channels x axis", "Engine + chamber wall diameter", "Channels width",
         "Channels height", "Channels area", "Gas viscosity", "Cp gas",
         "Gas conductivity", "Prandtl gaz", "Coeff Hg", "Sigma", " Twg ", " Twl ",
         "Heat flux", "Tl", "Reynolds CH4", "Coeff Hl", "Rho coolant",
         "Viscosity CH4", "Conductivity CH4", "Cp CH4", "Coolant velocity",
         "Coolant pressure", "Wall conductivity", "x real height", "y real height"))
    geometry1_writer.writerow(("x real height", "y real height"))
    geometry2_writer.writerow(("Engine + chamber wall radius", "x real height"))

    for i in range(0, nb_points):
        if i < nb_points_channel:
            geometry1_writer.writerow((newxhtre[i] * (-1000), newyhtre[i] * 1000))
            geometry2_writer.writerow((ycanaux[i] * 1000, newxhtre[i] * (-1000)))
            valuexport_writer.writerow((z_coord_list[i], r_coord_list[i], aire_saved[i], gamma_saved[i],
                                        mach_list_saved[i], static_pressure_list[i],
                                        hotgas_temperature_saved[i], z_coord_list[i], ycanaux[i],
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
                (z_coord_list[i], r_coord_list[i], aire_saved[i], gamma_saved[i],
                 mach_list_saved[i], static_pressure_list[i],
                 hotgas_temperature_saved[i], ' ', ' ', ' ', ' ', ' ', ' ',
                 ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                 ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))

    valuexport.close()
    geometry1.close()
    geometry2.close()

# %% Execution time display

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
print(f"█ Execution time for the initialisation       : {time_elapsed_i}                   █")
print("█                                                                          █")
print(f"█ Execution time for the main computation     : {time_elapsed_m}                   █")

if plot_detail >= 1:
    print("█                                                                          █")
    print(f"█ Execution time for the display of 1D        : {time_elapsed_d1}                   █")

if show_2D_temperature:
    print("█                                                                          █")
    print(f"█ Execution time for the display of 2D        : {time_elapsed_d2}                   █")

if do_final_3d_plot:
    print("█                                                                          █")
    print(f"█ Execution time for the display of 3D        : {time_elapsed_d3}                   █")

print("█                                                                          █")
print(f"█ Execution time for the end of the program   : {time_elapsed_e}                   █")
print("█                                                                          █")
print(f"█ Total execution time                        : {time_elapsed_t}                   █")
print("█                                                                          █")
print("███████████████████████████████████ END ████████████████████████████████████")

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
# Data
from Canaux import canaux

# Graphics
from heatequationsolve import carto2D
from volume3d import carto3d, view3d
import matplotlib.pyplot as plt
from tqdm import tqdm  # For progress bars

# Chemical species data
import fluid_properties as flp

start_time = time.perf_counter()  # Beginning of the timer

print("██████████████████████████ Cool The Engine V 2.0.0 █████████████████████████")
print("█                                                                          █")
print("█                  Innovative Propulsion Laboratory - IPL                  █")
print("█__________________________________________________________________________█")
print("█                                                                          █")
print("█ Initialisation                                                           █")

# %% Initial definitions

mesh_size = 0.25  # Distance between two points of calculation
x_coords_filename = f"input/{mesh_size}/x.txt"  # X coordinates of the Viserion
y_coords_filename = f"input/{mesh_size}/y.txt"  # Y coordinates of the Viserion
input_CEA_data = "input/Viserion_2023.txt"  # Viserion's parameters (found with CEA)

# Constant input_data_list
size2 = 16  # Used for the height of the display in 3D view
limitation = 0.05  # used to build the scales in 3D view
figure_dpi = 150  # Dots Per Inch (DPI) for all figures (lower=faster)
plot_detail = 1  # 0=No plots; 1=Important plots; 3=All plots
show_3d_plots = False
show_2D_temperature = False
do_final_3d_plot = False

# %% Reading input data
input_data_reader = csv.reader(open(input_CEA_data, "r"))
input_data_list = [row[1] for row in input_data_reader]

# Store CEA output in lists
sound_speed_init = float(input_data_list[0])  # Sound velocity in the chamber
sound_speed_throat = float(input_data_list[1])  # Sound velocity in the throat
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

# Store input dimensions in lists
curv_radius_pre_throat = float(input_data_list[12])  # Radius of curvature before the throat
curv_radius_after_throat = float(input_data_list[13])  # Radius of curvature after the throat
area_throat = float(input_data_list[14])  # Area at the throat
diam_throat = float(input_data_list[15])  # Throat diameter

# %% Import of the (X,Y) coordinates of the Viserion
x_coords_reader = csv.reader(open(x_coords_filename, "r"))
y_coords_reader = csv.reader(open(y_coords_filename, "r"))

# Storing the X,Y coordinates in lists
x_coord_list = [float(row[0]) / 1000 for row in x_coords_reader]
y_coord_list = [float(row[0]) / 1000 for row in y_coords_reader]
nb_points = len(x_coord_list)  # Number of points (or the index of the end of the divergent)

# Plot of the profile of the engine
if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(x_coord_list, y_coord_list, color='black')
    plt.title('Profile of the Viserion (left : chamber and right : divergent)', color='black')
    plt.show()

# Computation and plot of the mesh density of the engine
if plot_detail >= 3 and show_3d_plots:
    dist_between_pts = [abs(x_coord_list[i] - x_coord_list[i + 1]) for i in range(0, len(x_coord_list) - 1)]
    dist_between_pts.append(dist_between_pts[-1])
    colormap = plt.cm.binary
    inv = 1, 1, 1  # 1 means should be reversed
    view3d(inv, x_coord_list, y_coord_list, dist_between_pts, colormap, 'Mesh density (in m)', size2, limitation)

# %% Computation of the cross-sectional area along the engine
cross_section_area_list = [np.pi * r ** 2 for r in y_coord_list]

# Plots of the cross-sectionnal areas
if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(x_coord_list, cross_section_area_list, color='black')
    plt.title("Cross section area of engine (in m²) as a function of engine axis")
    plt.show()

# %% Adiabatic constant (gamma) parametrization
print("█                                                                          █")
print("█ Computing gamma                                                          █")
print("█                                                                          █")
i_conv = 0  # Index of the beginning of the convergent
y1 = 1
y2 = 1
while y1 == y2:  # Read y values two per two in order to detect the beginning of the convergent
    y1 = y_coord_list[i_conv]
    i_conv += 1
    y2 = y_coord_list[i_conv]

# Gamma in the cylindrical chamber
gamma_list = [gamma_c_input for i in range(0, i_conv)]  # Gamma is constant before the beginning of the convergent

# Gamma in the convergent
i_throat = y_coord_list.index(min(y_coord_list))  # Throat index
nb_points_convergent = i_throat - i_conv  # Number of points in the convergent
gamma_convergent = gamma_c_input
for m in range(-1, nb_points_convergent - 1):
    # Linear interpolation between beginning and end of convergent:
    # (yi+1)=((y2-y1)/(x2-x1))*abs((xi+1)-(xi))
    gamma_convergent += ((gamma_t_input - gamma_c_input) / (x_coord_list[i_throat] - x_coord_list[i_conv])) * abs(
        x_coord_list[i_conv + 1 + m] - x_coord_list[i_conv + m])
    gamma_list.append(gamma_convergent)

# Gamma in the divergent nozzle
nb_points_divergent = len(x_coord_list) - i_throat  # Number of points in the divergent
gamma_divergent = gamma_t_input
for q in range(-1, nb_points_divergent - 1):  # Linear interpolation between beginning and end of divergent
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

v_init_gas = (debit_LOX + debit_mass_coolant) / (rho_init * cross_section_area_list[0])  # Initial velocity of the gases
mach_init_gas = v_init_gas / sound_speed_init  # Initial mach number
mach_gas = mach_init_gas
mach_list = [mach_init_gas]

# Mach number computations along the engine
with tqdm(total=nb_points - 1,
          desc="█ Computing mach number        ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points - 1):
        mach_gas = t.mach_solv(cross_section_area_list[i], cross_section_area_list[i + 1],
                               mach_gas, gamma_list[i])
        mach_list.append(mach_gas)
        progressbar.update(1)

# Plots of the Mach number in the engine (2D/3D)
if plot_detail >= 1:
    plt.figure(dpi=figure_dpi)
    plt.plot(x_coord_list, mach_list, color='gold')
    plt.title("Mach number as a function of the engine axis")
    plt.show()

if plot_detail >= 1 and show_3d_plots:
    colormap = plt.cm.Spectral
    inv = 1, 1, 1  # 1 means should be reversed
    print("█ Plotting 3D graph                                                        █")
    print("█                                                                          █")
    view3d(inv, x_coord_list, y_coord_list, mach_list, colormap, 'Mach number of hot gases', size2, limitation)

# %% Static pressure computation
pressure_list = [Pc]  # (in Pa)

with tqdm(total=nb_points - 1,
          desc="█ Computing static pressure    ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points - 1):
        if i == nb_points + 1:  # If last point
            mach_gas = mach_list[i]
            mach_gas_2 = mach_list[i]
        else:  # All the other points
            mach_gas = mach_list[i]
            mach_gas_2 = mach_list[i + 1]
        pressure = t.pressure_solv(mach_gas, mach_gas_2, pressure_list[i], gamma_list[i])
        pressure_list.append(pressure)
        progressbar.update(1)

# Plot of the static pressure (2D/3D)
if plot_detail >= 2:
    plt.figure(dpi=figure_dpi)
    plt.plot(x_coord_list, pressure_list, color='gold')
    plt.title("Static pressure (in Pa) as a function of the engine axis")
    plt.show()

if plot_detail >= 2 and show_3d_plots:
    colormap = plt.cm.gist_rainbow_r
    inv = 1, 1, 1  # 1 means should be reversed
    print("█ Plotting 3D graph                                                        █")
    print("█                                                                          █")
    view3d(inv, x_coord_list, y_coord_list, pressure_list, colormap, 'Static pressure (in Pa)', size2, limitation)

# %% Hot gas temperature computation
hotgas_temp_list = [Tc]
with tqdm(total=nb_points - 1,
          desc="█ Computing gas temperature    ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points - 1):
        if i == nb_points + 1:
            mach_gas = mach_list[i]
            mach_gas_2 = mach_list[i]
        else:
            mach_gas = mach_list[i]
            mach_gas_2 = mach_list[i + 1]
        temperature = t.temperature_solv(mach_gas, mach_gas_2, hotgas_temp_list[i], gamma_list[i])
        hotgas_temp_list.append(temperature)
        progressbar.update(1)

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
    print("█ Plotting 3D graph                                                        █")
    view3d(inv, x_coord_list, y_coord_list, hotgas_temp_list, colormap, 'Temperature of the gases (in K)', size2,
           limitation)

# %% Dimensions
nbc = 40  # Number of channels
manifold_pos = 0.104  # Position of the manifold from the throat (in m)

# Widths
lrg_inj = 0.0045  # Width of the channel in at the injection plate (in m)
lrg_conv = 0.0025  # Width of the channel at the end of the cylindrical chamber (in m)
lrg_col = 0.0015  # Width of the channel in the throat (in m)
lrg_tore = 0.002  # Width of the channel at the manifold (in m)

# Heights
ht_inj = 0.002  # Height of the channel at the injection plate (in m)
ht_conv = 0.002  # Height of the channel at the end of the cylindrical chamber (in m)
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
if material == 0:  # Pure copper
    material_name = "pure copper"
elif material == 1:  # CuCrZr
    material_name = "cucrzr"
elif material == 2:
    material_name = "inconel"

# %% Properties of the coolant
fluid = "Methane"
density_cool_init = 425  # Density of the CH4 (kg/m^3)
Temp_cool_init = 110  # Initial temperature of the coolant (K)
debit_volum_coolant = debit_mass_coolant / density_cool_init  # Total volumic flow rate of the coolant (m^3/s)
Pressure_cool_init = 7000000  # Pressure of the coolant at inlet (Pa)
roughness = 15e-6  # Roughness (m)

# %% Computation of channel geometry
print("█                                                                          █")
print("█ Channel geometric computation                                            █")
print("█                                                                          █")
# Method 1
# xcanauxre, ycanauxre, larg_canalre, Areare, htre = canauxangl(x_coords_filename, y_coords_filename,
#                                                               nbc, lrg_col, ht_col, ht_c, ht_div, tore,
#                                                               debit_total, epaisseur_chemise)
# Method 2
# Pack the data in tuples
profile = (x_coord_list, y_coord_list)
widths = (lrg_inj, lrg_conv, lrg_col, lrg_tore)
heights = (ht_inj, ht_conv, ht_col, ht_tore)
thicknesses = (e_conv, e_col, e_tore)
coeffs = (n1, n2, n3, n4, n5, n6)

# Compute dimensions
xcanaux, ycanaux, larg_canal, larg_ailette_list, ht_canal, wall_thickness, area_channel, nb_points_channel \
    = canaux(profile, widths, heights, thicknesses, coeffs, manifold_pos, debit_volum_coolant, nbc, plot_detail,
             figure_dpi)

# Write the dimensions of the channels in a CSV file
file_name = "output/channelvalue.csv"
with open(file_name, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(("Engine x", "Engine y", "Channel width", "Rib width",
                     "Channel height", "Chamber wall thickness", "Channel area"))
    for i in range(0, nb_points_channel):
        writer.writerow((xcanaux[i], ycanaux[i], larg_canal[i], larg_ailette_list[i],
                         ht_canal[i], wall_thickness[i], area_channel[i]))

end_init_time = time.perf_counter()  # End of the initialisation timer
time_elapsed_i = time.ctime(end_init_time - start_time)[14:19]  # Initialisation time converted in min:sec
start_main_time = time.perf_counter()  # Start of the main solution timer

# %% Prepare the lists before main computation

wall_thickness.reverse()
xcanaux.reverse()
larg_canal.reverse()
area_channel.reverse()
ht_canal.reverse()
ycanaux.reverse()
# We reverse the lists in order to calculate from the manifold to the injection

# Save the data for exporting, before altering the original lists
hotgas_temperature_saved = hotgas_temp_list[:]
aire_saved = cross_section_area_list[:]
mach_list_saved = mach_list[:]
gamma_saved = gamma_list[:]

# Remove the data points before the manifold
hotgas_temp_list = hotgas_temp_list[:nb_points_channel]
cross_section_area_list = cross_section_area_list[:nb_points_channel]
mach_list = mach_list[:nb_points_channel]
gamma_list = gamma_list[:nb_points_channel]

gamma_list.reverse()
mach_list.reverse()
cross_section_area_list.reverse()
hotgas_temp_list.reverse()


# %% Main computation

def mainsolver(init_coolant_temp, init_coolant_pressure):
    """
    This is the main function used for solving the 1D case.
    The geometry is discretised into a 1 dimensionnal set of points.
    The function uses a marching algorithm, and computes all the relevant physical
    quantities at each point. The values obtained are then used on the next point.
    """

    # Lists containing the physical quantities at each point
    coolant_temp_list = [init_coolant_temp]
    coolant_pressure_list = [init_coolant_pressure]
    coolant_viscosity_list = [flp.viscosity(init_coolant_pressure, init_coolant_temp, fluid)]
    coolant_cond_list = [flp.conductivity(init_coolant_pressure, init_coolant_temp, fluid)]
    coolant_cp_list = [flp.cp(init_coolant_pressure, init_coolant_temp, fluid)]
    coolant_density_list = [flp.density(init_coolant_pressure, init_coolant_temp, fluid)]
    hotgas_viscosity_list = []
    hotgas_cp_list = []
    hotgas_cond_list = []
    hotgas_prandtl_list = []
    coolant_reynolds_list = []
    hg_list = []
    sigma_list = []
    wall_cond_list = []
    hotwall_temp_list = []
    coldwall_temp_list = []
    flux_list = []
    coolant_velocity_list = []
    sound_speed_list = []
    hl_normal_list = []
    hl_corrected_list = []
    hl_corrected_list_2 = []
    index_throat = ycanaux.index(min(ycanaux))

    length_from_inlet = 0.02
    coldwall_temp = 300
    hotwall_temp = 300

    with tqdm(total=nb_points_channel,
              desc=f"█ Global resolution            ",
              unit="|   █", bar_format="{l_bar}{bar}{unit}",
              ncols=76) as pbar_main:

        # Main computation loop
        for i in range(0, nb_points_channel):
            # Hydraulic diameter (4 * Area/Perimeter )
            Dhy = (2 * ht_canal[i] * larg_canal[i]) / (ht_canal[i] + larg_canal[i])

            # Velocity of the coolant
            v_cool = debit_mass_coolant / (nbc * coolant_density_list[i] * area_channel[i])
            coolant_velocity_list.append(v_cool)

            # Reynolds number of the coolant
            Re_cool = (v_cool * Dhy * coolant_density_list[i]) / coolant_viscosity_list[i]
            coolant_reynolds_list.append(Re_cool)

            # Prandtl number of the coolant
            Pr_cool = (coolant_viscosity_list[i] * coolant_cp_list[i]) / coolant_cond_list[i]

            # Compute viscosity, Cp, conductivity and Prandtl number of the hot gases
            hotgas_visc, hotgas_cp, hotgas_cond, hotgas_prandtl = t.hotgas_properties(hotgas_temp_list[i],
                                                                                      molar_mass,
                                                                                      gamma_list[i])

            # If last point in the list
            if i == len(xcanaux) - 1:
                # Distance between current point and the previous (Pythagoras)
                dl = ((xcanaux[i - 1] - xcanaux[i]) ** 2 + (ycanaux[i - 1] - ycanaux[i]) ** 2) ** 0.5
            else:
                # Distance between current point and the next (Pythagoras)
                dl = ((xcanaux[i + 1] - xcanaux[i]) ** 2 + (ycanaux[i + 1] - ycanaux[i]) ** 2) ** 0.5

            length_from_inlet += dl

            wall_cond = 350
            sigma = 1
            new_coldwall_temp = coldwall_temp + 10
            new_hotwall_temp = hotwall_temp + 10
            while abs(new_coldwall_temp - coldwall_temp) > 0.1 and abs(new_hotwall_temp - hotwall_temp) > 0.1:
                coldwall_temp = new_coldwall_temp
                hotwall_temp = new_hotwall_temp

                # Gas-side convective heat transfer coefficient (Bartz equation)
                hg = (0.0195 / (diam_throat ** 0.2) * (((hotgas_visc ** 0.2) * hotgas_cp) / (hotgas_prandtl ** 0.6)) * (
                        (Pc / c_star) ** 0.8) * ((diam_throat / curv_radius_pre_throat) ** 0.1) * (
                              (area_throat / cross_section_area_list[i]) ** 0.9)) * sigma

                # Coolant-side convective heat transfer coefficient from Taylor (NASA TN D-4332)
                Nu = 0.023 * Re_cool ** 0.705 * Pr_cool ** 0.8 * (coldwall_temp / coolant_temp_list[i]) ** -(
                        -0.57 - 1.59 * Dhy / (length_from_inlet))

                # Correction factor for the channel roughness
                xi = t.darcy_weisbach(Dhy, Re_cool, roughness) / t.darcy_weisbach(Dhy, Re_cool, 0)
                roughness_correction = xi * ((1 + 1.5 * Pr_cool ** (-1 / 6) * Re_cool ** (-1 / 8) * (Pr_cool - 1)) / (
                        1 + 1.5 * Pr_cool ** (-1 / 6) * Re_cool ** (-1 / 8) * (Pr_cool * xi - 1)))

                # Compute coolant-side convective heat-transfer coefficient
                hl = Nu * roughness_correction * (coolant_cond_list[i] / Dhy)

                # Fin dimensions
                D = 2 * (ycanaux[i] - wall_thickness[i])  # Diameter inside the engine
                fin_width = (np.pi * (D + ht_canal[i] + wall_thickness[i]) - nbc * larg_canal[
                    i]) / nbc  # Width of the fin

                # Correct for the fin effect (unknown source)
                m_ = ((2 * hl) / (fin_width * wall_cond)) ** 0.5
                hl_cor = hl * ((nbc * larg_canal[i]) / (np.pi * D)) + nbc * (
                        (2 * hl * wall_cond * (((np.pi * D) / nbc) - larg_canal[i])) ** 0.5) * (
                                 (np.tanh(m_ * ht_canal[i])) / (np.pi * D))

                # Correct for the fin effect (Luka Denies)
                intermediate_calc_1 = ((2 * hl * fin_width) / wall_cond) ** 0.5 * ht_canal[i] / fin_width
                nf = np.tanh(intermediate_calc_1) / intermediate_calc_1
                hl_cor2 = hl * (larg_canal[i] + 2 * nf * ht_canal[i]) / (larg_canal[i] + fin_width)

                # Alternative but equivalent way of computing the heat flux and wall temperatures (Luka Denies)
                flux = (hotgas_temp_list[i] - coolant_temp_list[i]) / (
                        1 / hg + 1 / hl_cor + wall_thickness[i] / wall_cond)
                new_hotwall_temp = hotgas_temp_list[i] - flux / hg
                new_coldwall_temp = coolant_temp_list[i] + flux / hl

                # Compute sigma (used in the Bartz equation)
                T_hotgas_throat = hotgas_temp_list[index_throat]
                mach_hot_gases = mach_list[i]
                sigma = (((new_hotwall_temp / (2 * T_hotgas_throat)) * (
                        1 + (((gamma_list[i] - 1) / 2) * (mach_hot_gases ** 2))) + 0.5) ** -0.68) * (
                                (1 + (((gamma_list[i] - 1) / 2) * (mach_hot_gases ** 2))) ** -0.12)

                # Compute thermal conductivity of the solid at a given temperature
                wall_cond = t.conductivity(Twg=new_hotwall_temp, Twl=new_coldwall_temp, material_name=material_name)

            coldwall_temp = new_coldwall_temp
            hotwall_temp = new_hotwall_temp

            # Compute heat exchange area between two points
            # Cold-wall version (Julien)
            dA_1 = 2 * dl * (larg_canal[i] + ht_canal[i])
            # Hot-wall version (Luka Denies)
            dA_2 = (np.pi * D * dl) / nbc

            # New temperature at next point
            delta_T_coolant = ((flux * dA_1) / ((debit_mass_coolant / nbc) * coolant_cp_list[i]))
            new_coolant_temp = coolant_temp_list[i] + delta_T_coolant

            # Solving Colebrook's formula to obtain the Darcy-Weisbach friction factor
            frict_factor = t.darcy_weisbach(Dhy, Re_cool, roughness)

            # Computing pressure loss with the Darcy-Weisbach friction factor
            delta_p = 0.5 * frict_factor * (dl / Dhy) * coolant_density_list[i] * v_cool ** 2
            new_coolant_pressure = coolant_pressure_list[i] - delta_p

            # Computing the new properties of the CH4
            if new_coolant_pressure < 0:
                raise ValueError("Negative pressure ! Pressure drop is too high.")
            new_cool_visc = flp.viscosity(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)
            new_cool_cond = flp.conductivity(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)
            new_cool_cp = flp.cp(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)
            new_cool_dens = flp.density(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)
            new_cool_sound_spd = flp.sound_speed(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)

            # Store the results
            hotgas_viscosity_list.append(hotgas_visc)
            hotgas_cp_list.append(hotgas_cp)
            hotgas_cond_list.append(hotgas_cond)
            hotgas_prandtl_list.append(hotgas_prandtl)
            hg_list.append(hg)
            hl_normal_list.append(hl)
            hl_corrected_list.append(hl_cor)
            hl_corrected_list_2.append(hl_cor2)
            hotwall_temp_list.append(coldwall_temp)
            coldwall_temp_list.append(hotwall_temp)
            flux_list.append(flux)
            sigma_list.append(sigma)
            wall_cond_list.append(wall_cond)
            coolant_pressure_list.append(new_coolant_pressure)
            coolant_temp_list.append(new_coolant_temp)
            coolant_viscosity_list.append(new_cool_visc)
            coolant_cond_list.append(new_cool_cond)
            coolant_cp_list.append(new_cool_cp)
            coolant_density_list.append(new_cool_dens)
            sound_speed_list.append(new_cool_sound_spd)

            pbar_main.update(1)

        return hl_corrected_list, hl_corrected_list_2, hotgas_viscosity_list, \
               hotgas_cp_list, hotgas_cond_list, hotgas_prandtl_list, hg_list, \
               hotwall_temp_list, coldwall_temp_list, flux_list, sigma_list, \
               coolant_reynolds_list, coolant_temp_list, coolant_viscosity_list, \
               coolant_cond_list, coolant_cp_list, coolant_density_list, \
               coolant_velocity_list, coolant_pressure_list, wall_cond_list, \
               sound_speed_list, hl_normal_list


# First iteration of the solving
hlcor_list, hlcor_list_2, hotgas_visc_list, hotgas_cp_list, hotgas_cond_list, hotgas_prandtl_list, hg_list, \
hotwall_temp_list, coldwall_temp_list, flux_list, \
sigma_list, coolant_reynolds_list, tempcoolant_list, visccoolant_list, \
condcoolant_list, cpcoolant_list, densitycoolant_list, velocitycoolant_list, \
pcoolant_list, wallcond_list, sound_speed_coolant_list, hlnormal_list \
    = mainsolver(Temp_cool_init, Pressure_cool_init)

end_m = time.perf_counter()  # End of the main solution timer
time_elapsed_m = time.ctime(end_m - start_main_time)[14:19]  # Main elapsed time converted in minutes:secondes
start_d2 = time.perf_counter()  # Start of the display of 2D timer

# %% Display of the 1D analysis results
mach_03 = [x * 0.3 for x in sound_speed_coolant_list]

if plot_detail >= 1:
    plt.figure(dpi=figure_dpi)
    plt.plot(xcanaux, hlcor_list_2, color='blue', label='Hl corrected (Luka Denies)')
    plt.plot(xcanaux, hlcor_list, color='blue', label='Hl corrected (Julien)')
    plt.plot(xcanaux, hlnormal_list, color='cyan', label='Hl')
    plt.title("Convection coeff as a function of the engine axis")
    plt.legend(loc='upper left')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(xcanaux, coldwall_temp_list, color='red', label='Twg')
    plt.plot(xcanaux, hotwall_temp_list, color='blue', label='Twl')
    plt.title('Wall temperature (in K) as a function of engine axis')
    plt.legend(loc='lower left')
    plt.show()

    tempcoolant_list.pop()
    plt.figure(dpi=figure_dpi)
    plt.plot(xcanaux, tempcoolant_list, color='blue')
    plt.title('Coolant temperature (in K) as a function of engine axis')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(xcanaux, velocitycoolant_list, color='blue', label='Coolant')
    plt.plot(xcanaux, mach_03, color='orange', label='Mach 0.3 limit')
    plt.title('Velocity (in m/s) of the coolant as a function of engine axis')
    plt.legend(loc='upper left')
    plt.show()

    pcoolant_list.pop()
    plt.figure(dpi=figure_dpi)
    plt.plot(xcanaux, pcoolant_list, color='orange')
    plt.title('Pressure drop in the cooling channels')
    plt.show()

if plot_detail >= 2:
    wallcond_list.pop()
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

if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(xcanaux, coolant_reynolds_list, color='blue')
    plt.title("Reynolds number of the coolant as a function of the engine axis")
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
    plt.title('Sound velocity of the coolant (in m/s) as a function of engine axis')
    plt.show()

if plot_detail >= 1 and show_3d_plots:
    colormap = plt.cm.plasma
    inv = 0, 0, 0
    view3d(inv, xcanaux, ycanaux, flux_list, colormap, "Heat flux (in MW/m²)", size2, limitation)

    colormap = plt.cm.coolwarm
    inv = 0, 0, 0
    view3d(inv, xcanaux, ycanaux, tempcoolant_list, colormap, "Temperature of the coolant (in K)", size2, limitation)

if plot_detail >= 2 and show_3d_plots:
    colormap = plt.cm.magma
    inv = 0, 0, 0  # 1 means should be reversed
    view3d(inv, xcanaux, ycanaux, coldwall_temp_list, colormap, "Wall temperature on the gas side (in K)", size2,
           limitation)

# %% Flux computation in 2D and 3D
print("█                                                                          █")
print("█ Display of results                                                       █")
print("█                                                                          █")
"""2D flux computation"""
larg_ailette_list.reverse()

if show_2D_temperature:
    # At the beginning of the chamber
    print("█ Results at the beginning of the chamber :                                █")
    pas = larg_ailette_list[-1] + larg_canal[-1]
    epaisseur = e_conv
    hauteur = ht_canal[-1]
    largeur = larg_canal[-1]
    Hg = hg_list[-1]
    Tg = hotgas_temp_list[-1]
    Hl = hlnormal_list[-1]
    Tl = tempcoolant_list[-1]
    dx = 0.00004  # *3.5
    wall_cond_throat = wallcond_list[-1]
    where = " at the beginning of the chamber"
    t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, wall_cond_throat, Tg, Hl, Tl, 5, 1, 1, where,
                  show_2D_temperature)

    # At the throat
    print("█ Results at the throat :                                                  █")
    pos_col = ycanaux.index(min(ycanaux))
    pas = larg_ailette_list[pos_col] + larg_canal[pos_col]
    epaisseur = e_col
    hauteur = ht_canal[pos_col]
    largeur = larg_canal[pos_col]
    Hg = hg_list[pos_col]
    Tg = hotgas_temp_list[pos_col]
    Hl = hlnormal_list[pos_col]
    Tl = tempcoolant_list[pos_col]
    dx = 0.000025  # *3.5
    wall_cond_throat = wallcond_list[pos_col]
    where = " at the throat"
    t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, wall_cond_throat, Tg, Hl, Tl, 15, 1, 2, where,
                  show_2D_temperature)

    # At the end of the divergent
    print("█ Results at the end of the divergent :                                    █")
    pas = larg_ailette_list[0] + larg_canal[0]
    epaisseur = e_tore
    hauteur = ht_canal[0]
    largeur = larg_canal[0]
    Hg = hg_list[0]
    Tg = hotgas_temp_list[0]
    Hl = hlnormal_list[0]
    Tl = tempcoolant_list[0]
    dx = 0.00004
    wall_cond_throat = wallcond_list[0]
    where = " at the end of the divergent"
    t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, wall_cond_throat, Tg, Hl, Tl, 5, 1, 1, where,
                  show_2D_temperature)

end_d2 = time.perf_counter()  # End of the display of 2D timer
time_elapsed_d2 = time.ctime(end_d2 - start_d2)[14:19]  # Display of 2D elapsed time converted in minutes:secondes

"Computation for 3D graph"
if do_final_3d_plot:
    # 3D display
    start_3d = time.perf_counter()
    eachT = []
    lim1 = 0
    lim2 = 650
    dx = 0.0001
    with tqdm(total=nb_points_channel,
              desc="█ 3D graph initialisation      ",
              unit="|   █", bar_format="{l_bar}{bar}{unit}",
              ncols=76) as progressbar:
        for i in range(0, nb_points_channel, 1):
            # if lim1 <= i <= lim2:
            #     dx = 0.0001
            # else:
            #     dx = 0.0001

            wall_cond_throat = wallcond_list[i]
            t3d = carto2D(larg_ailette_list[i] + larg_canal[i], wall_thickness[i], ht_canal[i], larg_canal[i], dx,
                          hg_list[i], wall_cond_throat,
                          hotgas_temp_list[i], hlnormal_list[i], tempcoolant_list[i], 3, 0, 1, "", plot_detail)
            eachT.append(t3d)
            progressbar.update(1)

    carto3d([0, 0, 0], xcanaux, ycanaux, eachT, plt.cm.Spectral_r, '3D view of wall temperatures (in K)', nbc,
            limitation)

    end_d3 = time.perf_counter()  # End of the display of 3D timer
    time_elapsed_d3 = time.ctime(end_d3 - start_d2)[14:19]  # Display of 3D elapsed time converted in minutes:secondes
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

# %% Preparation of the lists for CAD modelisation
"Changing the coordinates of the height of the channels (otherwise it is geometrically wrong)"

angles = []
newxhtre = []
newyhtre = []
print("█ Computing channel height                                                 █")
for i in range(0, nb_points_channel):
    if i == 0:
        angle = 0
        angles.append(angle)
    elif i == (nb_points_channel - 1):
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
print("█ Checking channel height                                                  █")
for i in range(0, nb_points_channel):
    verifhtre = (((newxhtre[i] - xcanaux[i]) ** 2) + ((newyhtre[i] - ycanaux[i]) ** 2)) ** 0.5
    verification.append(verifhtre)

if plot_detail >= 3:
    plt.figure(dpi=figure_dpi)
    plt.plot(newxhtre, newyhtre, color='blue', label='New height')
    plt.plot(xcanaux, ycanaux, color='chocolate', label='Former height')
    plt.title("Geometrical aspect of the channel (height as a function of the engine axis)")
    plt.axis("equal")
    plt.legend(loc='upper left')
    plt.show()

    plt.figure(dpi=figure_dpi)
    plt.plot(xcanaux, verification)
    plt.title("Checking the height of the generated channels")
    plt.show()

# %% Writing the results of the study in a CSV file
print("█ Writing results in .csv files                                            █")

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
geometry2_writer.writerow(("Engine + chamber wall radius", "x real height"))

for i in range(0, nb_points):
    if i < nb_points_channel:
        geometry1_writer.writerow((newxhtre[i] * (-1000), newyhtre[i] * 1000))
        geometry2_writer.writerow((ycanaux[i] * 1000, newxhtre[i] * (-1000)))
        valuexport_writer.writerow((x_coord_list[i], y_coord_list[i], aire_saved[i], gamma_saved[i],
                                    mach_list_saved[i], pressure_list[i],
                                    hotgas_temperature_saved[i], xcanaux[i], ycanaux[i],
                                    larg_canal[i], ht_canal[i], area_channel[i], hotgas_visc_list[i],
                                    hotgas_cp_list[i], hotgas_cond_list[i], hotgas_prandtl_list[i],
                                    hg_list[i], sigma_list[i], coldwall_temp_list[i],
                                    hotwall_temp_list[i], flux_list[i], tempcoolant_list[i],
                                    coolant_reynolds_list[i], hlnormal_list[i], densitycoolant_list[i],
                                    visccoolant_list[i],
                                    condcoolant_list[i], cpcoolant_list[i], velocitycoolant_list[i], pcoolant_list[i],
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
time_elapsed_e = time.ctime(end_t - start_e)[14:19]  # End elapsed time converted in minutes:secondes
time_elapsed_t = time.ctime((end_t - start_e) + (end_d2 - start_time))[14:19]
# Total elapsed time (without time waited to choose 3D) converted in minutes:secondes

print("█                                                                          █")
print("█__________________________________________________________________________█")
print("█                                                                          █")
print("█ Execution time for the initialisation       :", time_elapsed_i, "                     █")
print("█                                                                          █")
print("█ Execution time for the computation of main  :", time_elapsed_m, "                     █")
print("█                                                                          █")
print("█ Execution time for the display of 2D        :", time_elapsed_d2, "                     █")
print("█                                                                          █")
print("█ Execution time of the end of the program    :", time_elapsed_e, "                     █")

if do_final_3d_plot:
    time_elapsed_t_w3D = time.ctime((end_t - start_3d) + (end_d2 - start_time))[14:19]
    # Total elapsed time with 3D computation (without time waited to choose 3D) converted in minutes:secondes
    print("█                                                                          █")
    print("█ Execution time for the display of 3D        :", time_elapsed_d3, "                     █")
    print("█                                                                          █")
    print("█ Total execution time with 3D computation    :", time_elapsed_t_w3D, "                     █")

print("█                                                                          █")
print("█ Total execution time without 3D computation :", time_elapsed_t, "                     █")
print("█                                                                          █")
print("███████████████████████████████████ END ████████████████████████████████████")

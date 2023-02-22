"""
Created on Fri Nov 27 14:47:27 2020

Original author: Julien S

Refactored and improved by Mehdi D, Paul B, Paul M, Eve X and Antoine R
"""

import time
import csv

# Calculations
import numpy as np
from machsolve import mach_solv
from pressuresolve import pressure_solv
from temperaturesolve import temperature_solv
import scipy.optimize as opt

# Data
from Canaux import canaux

# Graphics
from heatequationsolve import carto2D
from volume3d import carto3d, view3d
import matplotlib.pyplot as plt
from tqdm import tqdm  # For progress bars

# Chemical species data
from CoolProp.CoolProp import PropsSI
import methan as meth


# %% Function definitions
def tempcorrige(temp_original, gamma, mach_number):
    Pr = 4 * gamma / (9 * gamma - 5)
    temp_corrected = temp_original * (1 + (Pr ** 0.33) * (
            (gamma - 1) / gamma) * mach_number ** 2) / (1 + (
            (gamma - 1) / gamma) * mach_number ** 2)
    return temp_corrected


def flux_equations(vars, *data):
    """
    Used by scipy.optimize.fsolve() to compute hot and cold wall temperature.
    """
    t_hot, t_cold = vars  # Initial guess
    hg, hl, t_g, t_c, wall_cond, wall_thickness = data

    # System of equations to solve
    f1 = hg * (t_g - t_hot) - (wall_cond / wall_thickness) * (t_hot - t_cold)
    f2 = hl * (t_cold - t_c) - (wall_cond / wall_thickness) * (t_hot - t_cold)
    return [f1, f2]


def conductivity(Twg: float, Twl: float, material_name: str):
    T_avg = (Twg + Twl) * 0.5
    if material_name == "pure copper":
        return -0.065665 * T_avg + 421.82
    if material_name == "cucrzr":
        return -0.0269 * T_avg + 365.33


def mu(gas_temp, molar_mass_, gamma):
    dyn_viscosity = 17.858 * (46.61 * 10 ** (-10)) * (molar_mass_ ** 0.5) * ((9 / 5) * gas_temp) ** 0.6
    Cp = 8314 * gamma / ((gamma - 1) * molar_mass_)
    Lamb = dyn_viscosity * (Cp + (8314.5 / (4 * molar_mass_)))
    Pr = 4 * gamma / (9 * gamma - 5)

    return dyn_viscosity, Cp, Lamb, Pr


start_time = time.perf_counter()  # Beginning of the timer

print("██████████████████████████ Cool The Engine V 2.0.0 █████████████████████████")
print("█                                                                          █")
print("█                  Innovative Propulsion Laboratory - IPL                  █")
print("█                                                                          █")
print("█ Initialisation                                                           █")

# %% Initial definitions

mesh_size = 0.25  # Distance between two points of calculation
x_coords_filename = f"input/{mesh_size}/x.txt"  # X coordinates of the Viserion
y_coords_filename = f"input/{mesh_size}/y.txt"  # Y coordinates of the Viserion
input_CEA_data = "input/Viserion_2023.txt"  # Viserion's parameters (found with CEA)

# Constant input_data_list
size2 = 16  # Used for the height of the display in 3D view
machtype = 0  # 0 for one equation and else for another equation in calculation of mach
limitation = 0.05  # used to build the scales in 3D view

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

# Plot of the profile of the engine
"""
plt.figure(dpi=300)
plt.plot(x_coord_list, y_coord_list, color='black')
plt.title('Profile of the Viserion (left : chamber and right : divergent', color='black')
plt.show()
"""

# Computation and plot of the mesh density of the engine
dist_between_pts = [abs(x_coord_list[i] - x_coord_list[i + 1]) for i in range(0, len(x_coord_list) - 1)]
dist_between_pts.append(dist_between_pts[-1])
colormap = plt.cm.binary
inv = 1, 1, 1  # 1 means should be reversed
view3d(inv, x_coord_list, y_coord_list, dist_between_pts, colormap, 'Mesh density', size2, limitation)

# %% Computation of the cross-sectional area along the engine
cross_section_area_list = [np.pi * r ** 2 for r in y_coord_list]

# Plots of the cross-sectionnal areas
"""
plt.figure(dpi=300)
plt.plot(x_value, aire, color='black')
plt.title("Area (in m²) as a function of engine axis")
plt.show()
"""
print("█                                                                          █")

# %% Adiabatic constant (gamma) parametrization
print("█ Computing gamma                                                          █")
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
plt.figure(dpi=300)
plt.plot(x_coord_list, gamma_list, color='gold')
plt.title("Gamma linearisation as a function of engine axis")
plt.show()

print("█                                                                          █")

# %% Mach number computation
"Computation of the initial velocity and mach number of the gases"
v_init_gas = (debit_LOX + debit_mass_coolant) / (rho_init * cross_section_area_list[0])  # Initial velocity of the gases
mach_init_gas = v_init_gas / sound_speed_init  # Initial mach number
mach_gas = mach_init_gas
mach_list = [mach_init_gas]
nb_points = len(x_coord_list)  # Number of points (or the index of the end of the divergent)

"Mach number computations along the engine"
with tqdm(total=nb_points - 1,
          desc="█ Computing mach number        ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points - 1):
        mach_gas = mach_solv(cross_section_area_list[i],
                             cross_section_area_list[i + 1],
                             mach_gas, gamma_list[i], i, machtype)
        mach_list.append(mach_gas)
        progressbar.update(1)

# Plots of the Mach number in the engine (2D/3D)
plt.figure(dpi=300)
plt.plot(x_coord_list, mach_list, color='gold')
plt.title("Mach number as a function of the engine axis")
plt.show()

colormap = plt.cm.Spectral
inv = 1, 1, 1  # 1 means should be reversed
print("█ Plotting 3D graph                                                        █")
print("█                                                                          █")
view3d(inv, x_coord_list, y_coord_list, mach_list, colormap, 'Mach number', size2, limitation)

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
        pressure = pressure_solv(mach_gas, mach_gas_2, pressure_list[i], gamma_list[i])
        pressure_list.append(pressure)
        progressbar.update(1)

# Plot of the static pressure (2D/3D)
plt.figure(dpi=300)
plt.plot(x_coord_list, pressure_list, color='gold')
plt.title("Static pressure as a function of the engine axis")
plt.show()

colormap = plt.cm.gist_rainbow_r
inv = 1, 1, 1  # 1 means should be reversed
print("█ Plotting graph                                                           █")
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
        temperature = temperature_solv(mach_gas, mach_gas_2, hotgas_temp_list[i], gamma_list[i])
        hotgas_temp_list.append(temperature)
        progressbar.update(1)

# List of corrected gas temperatures (max diff with original is about 75 K)
hotgas_temp_list_corrected = [tempcorrige(hotgas_temp_list[i], gamma_list[i], mach_list[i]) for i in
                              range(0, nb_points)]

# Plots of the temperature in the engine (2D/3D)
"""
plt.figure(dpi=300)
plt.plot(x_coord_list, hotgas_temp_list, color='gold')
plt.title("Gas temperature (in K) as a function of the engine axis")
plt.show()
colormap = plt.cm.coolwarm
inv = 1, 1, 1  # 1 means should be reversed
print("█ Plotting graph...                                                        █")
view3d(inv, x_coord_list, y_coord_list, hotgas_temp_list, colormap, 'Temperature of the gases (in K)', size2,
       limitation)
"""
# %% Dimensions
nbc = 40  # Number of channels
manifold_pos = 0.104  # Position of the manifold from the throat (in m)

# Widths
lrg_inj = 0.002  # Width of the channel in at the injection plate (in m)
lrg_conv = 0.002  # Width of the channel at the end of the cylindrical chamber (in m)
lrg_col = 0.002  # Width of the channel in the throat (in m)
lrg_tore = 0.002  # Width of the channel at the manifold (in m)

# Heights
ht_inj = 0.002  # Height of the channel at the injection plate (in m)
ht_conv = 0.002  # Height of the channel at the end of the cylindrical chamber (in m)
ht_col = 0.002  # Height of the channel in the throat (in m)
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
material = 0
if material == 0:  # Pure copper
    material_name = "pure copper"
elif material == 1:  # CuCrZr
    material_name = "cucrzr"

# %% Properties of the coolant
fluid = "Methane"
density_cool_init = 425  # Density of the CH4 (kg/m^3)
Temp_cool_init = 110  # Initial temperature of the coolant (K)
debit_volum_coolant = debit_mass_coolant / density_cool_init  # Total volumic flow rate of the coolant (m^3/s)
Pressure_cool_init = 7000000  # Pressure of the coolant at inlet (Pa)
roughness = 15e-6  # Roughness (m)

# %% Computation of channel geometry
print("█ Channel geometric computation                                            █")

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
xcanaux, ycanaux, larg_canal, larg_ailette, ht_canal, wall_thickness, area_channel, nb_points_channel \
    = canaux(profile, widths, heights, thicknesses, coeffs, manifold_pos, debit_volum_coolant, nbc)

# Write the dimensions of the channels in a CSV file
file_name = "channelvalue.csv"
with open(file_name, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(("Engine x", "Engine y", "Channel width", "Rib width",
                     "Channel height", "Chamber wall thickness", "Channel area"))
    for i in range(0, nb_points_channel):
        writer.writerow((xcanaux[i], ycanaux[i], larg_canal[i], larg_ailette[i],
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
def mainsolver(sigma_list, coolant_density_list, coolant_temp_list,
               coolant_viscosity_list, coolant_cond_list, coolant_cp_list,
               coolant_pressure_list, wall_cond_list, current_rep, total_reps):
    """
    This is the main function used for solving the 1D case.
    The geometry is discretised into a 1 dimensionnal set of points.
    The function uses a marching algorithm, and computes all the relevant physical
    quantities at each point. The values obtained are then used on the next point.
    """
    # Lists containing the physical quantities at each point
    hotgas_viscosity_list = []
    hotgas_cp_list = []
    coolant_reynolds_list = []
    hotgas_cond_list = []
    hotgas_prandtl_list = []
    hg_list = []
    hotwall_temperature_list = [110]
    coldwall_temperature_list = [110]
    flux_density_list = []
    coolant_velocity_list = []
    sound_speed_list = []
    hl_normal_list = []
    hl_corrected_list = []
    singpertes = [coolant_pressure_list[0]]
    Pcoolant2 = [coolant_pressure_list[0]]
    index_throat = ycanaux.index(min(ycanaux))

    length_from_inlet = 0
    with tqdm(total=nb_points_channel,
              desc=f"█ Global resolution {current_rep}/{total_reps}        ",
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
            hotgas_visc, hotgas_cp, hotgas_cond, hotgas_prandtl = mu(hotgas_temp_list[i], molar_mass, gamma_list[i])

            # Store the results
            hotgas_viscosity_list.append(hotgas_visc)
            hotgas_cp_list.append(hotgas_cp)
            hotgas_cond_list.append(hotgas_cond)
            hotgas_prandtl_list.append(hotgas_prandtl)

            # Gas-side convective heat transfer coefficient (Bartz equation)
            hg = (0.026 / (diam_throat ** 0.2) * (((hotgas_visc ** 0.2) * hotgas_cp) / (hotgas_prandtl ** 0.6)) * (
                    (Pc / c_star) ** 0.8) * ((diam_throat / curv_radius_pre_throat) ** 0.1) * (
                          (area_throat / cross_section_area_list[i]) ** 0.9)) * sigma_list[
                     i]
            hg_list.append(hg)

            # Radiative heat flux assuming black body radiation
            # Tg = hotgas_temperature[i]
            # steff = 5.6697 * 10 ** (-8)  # Stefan-Boltzmann constant
            # emissivity = 0.02  # 2% emissivity for CH4
            # qr = emissivity * steff * (T1 ** 4)  # Radiative heat flux

            # If last point in the list
            if i == len(xcanaux) - 1:
                # Distance between current point and the previous (Pythagoras)
                dr = ((xcanaux[i - 1] - xcanaux[i]) ** 2 + (ycanaux[i - 1] - ycanaux[i]) ** 2) ** 0.5
            else:
                # Distance between current point and the next (Pythagoras)
                dr = ((xcanaux[i + 1] - xcanaux[i]) ** 2 + (ycanaux[i + 1] - ycanaux[i]) ** 2) ** 0.5

            length_from_inlet += dr

            # Modified correlation from Taylor (NASA TN D-4332)
            inlet_buffer = 0.02
            Nu = 0.023 * Re_cool ** 0.705 * Pr_cool ** 0.8 * (
                    coldwall_temperature_list[i] / coolant_temp_list[i]) ** -(
                    -0.57 - 1.59 * Dhy / (length_from_inlet + inlet_buffer))

            # Modified Pizzarelli correlation (best model)
            # lambda_pizza = None  # This requires properties of the coolant in the bulk flow AND near the wall
            # Nu = 0.0082 * Re_hg ** 0.8 * hotgas_prandtl ** 0.16 * lambda_pizza ** 0.33 * ( 1 + 10 * (Dhy / length_from_inlet)) ** 0.5

            # Compute coolant-side convective heat-transfer coefficient
            hl = Nu * (coolant_cond_list[i] / Dhy)
            hl_normal_list.append(hl)

            # Fin dimensions
            D = 2 * (ycanaux[i] - wall_thickness[i])
            d = (np.pi * (D + ht_canal[i] + wall_thickness[i]) - nbc * larg_canal[i]) / nbc
            m = float(((2 * hl) / (d * wall_cond_list[i])) ** 0.5)

            # Correct for the fin effect
            hl_cor = hl * ((nbc * larg_canal[i]) / (np.pi * D)) + nbc * \
                     ((2 * hl * wall_cond_list[i] * (((np.pi * D) / nbc) - larg_canal[i])) ** 0.5) * (
                             (np.tanh(m * ht_canal[i])) / (np.pi * D))
            hl_corrected_list.append(hl_cor)

            # Solve a system to obtain the temperatures in the hot and cold wall
            T_hotwall, T_coldwall = opt.fsolve(func=flux_equations,
                                               x0=np.array([900.0, 700.0]),
                                               args=(hg, hl_cor,
                                                     hotgas_temp_list[i],
                                                     coolant_temp_list[i],
                                                     wall_cond_list[i],
                                                     wall_thickness[i]))

            # Temperature at the walls
            hotwall_temperature_list.append(T_hotwall)
            coldwall_temperature_list.append(T_coldwall)

            # Store heat flux density through the coolant side (W/m²)
            flux_density = hl * (coldwall_temperature_list[i] - coolant_temp_list[i]) * 0.000001
            flux_density_list.append(flux_density)

            # Compute sigma (used in the Bartz equation)
            Twg = hotwall_temperature_list[i]
            Twl = coldwall_temperature_list[i]
            T_hotgas_throat = hotgas_temp_list[index_throat]
            Mak = mach_list[i]
            sigma = 1 / ((((Twg / (2 * T_hotgas_throat)) * (
                    1 + (((gamma_list[i] - 1) / 2) * (Mak ** 2))) + 0.5) ** 0.68) * (
                                 (1 + (((gamma_list[i] - 1) / 2) * (Mak ** 2))) ** 0.12))
            sigma_list.append(sigma)

            # Compute thermal conductivity of the solid at a given temperature
            wall_cond = conductivity(Twg=Twg, Twl=Twl, material_name="cucrzr")
            wall_cond_list.append(wall_cond)

            # Compute heat exchange area between two points
            dA = dr * (2 * larg_canal[i] + 2 * ht_canal[i])

            # Compute total heat flux through the surface between 2 points
            Q = flux_density * dA * 1000000

            # Mass flow rate through one channel
            debit_mass_par_canal = debit_mass_coolant / nbc

            # New temperature at next point
            delta_T_coolant = (Q / (debit_mass_par_canal * coolant_cp_list[i]))
            coolant_temp_list.append(coolant_temp_list[i] + delta_T_coolant)

            # Solving Colebrook's formula to obtain the Darcy-Weisbach friction factor
            friction_factor_1 = 1e-3
            friction_factor_2 = (1 / (-2 * np.log10(
                ((roughness / (Dhy * 3.7)) + 2.51 / (Re_cool * (friction_factor_1 ** 0.5)))))) ** 2

            while abs((friction_factor_1 / friction_factor_2) - 1) > 0.0000001:
                friction_factor_1 = friction_factor_2
                friction_factor_2 = (1 / (-2 * np.log10(
                    ((roughness / (Dhy * 3.7)) + 2.51 / (Re_cool * (friction_factor_1 ** 0.5)))))) ** 2

            # Computing pressure loss with the Darcy-Weisbach friction factor
            delta_p = 0.5 * friction_factor_2 * (dr / Dhy) * coolant_density_list[i] * v_cool ** 2

            # Compute new pressure
            coolant_pressure_list.append(coolant_pressure_list[i] - delta_p)

            # Computing the new properties of the CH4
            new_cool_visc = PropsSI("V", "P", coolant_pressure_list[i + 1], "T", coolant_temp_list[i], fluid)
            new_cool_cond = PropsSI("L", "P", coolant_pressure_list[i + 1], "T", coolant_temp_list[i], fluid)
            new_cool_cp = PropsSI("C", "P", coolant_pressure_list[i + 1], "T", coolant_temp_list[i], fluid)
            new_cool_dens = PropsSI("D", "P", coolant_pressure_list[i + 1], "T", coolant_temp_list[i], fluid)
            new_cool_sound_spd = PropsSI("A", "P", coolant_pressure_list[i + 1], "T", coolant_temp_list[i], fluid)

            coolant_viscosity_list.append(new_cool_visc)
            coolant_cond_list.append(new_cool_cond)
            coolant_cp_list.append(new_cool_cp)
            coolant_density_list.append(new_cool_dens)
            sound_speed_list.append(new_cool_sound_spd)

            pbar_main.update(1)

        return hl_corrected_list, hotgas_viscosity_list, hotgas_cp_list, \
               hotgas_cond_list, hotgas_prandtl_list, hg_list, \
               hotwall_temperature_list, coldwall_temperature_list, \
               flux_density_list, sigma_list, coolant_reynolds_list, \
               coolant_temp_list, coolant_viscosity_list, coolant_cond_list, \
               coolant_cp_list, coolant_density_list, coolant_velocity_list, \
               coolant_pressure_list, wall_cond_list, sound_speed_list, \
               hl_normal_list, singpertes, Pcoolant2


Sig = [1]
LambdaTC = [350]
Tcoolant = [Temp_cool_init]
Pcoolant = [Pressure_cool_init]
visccoolant = [meth.viscCH4(Pressure_cool_init, Temp_cool_init, fluid)]
condcoolant = [meth.condCH4(Pressure_cool_init, Temp_cool_init, fluid)]
Cpmeth = [meth.cpCH4(Pressure_cool_init, Temp_cool_init, fluid)]
rho = [meth.densityCH4(Pressure_cool_init, Temp_cool_init, fluid)]
reps_max = 2

# First iteration of the solving
hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, \
inwall_temperature, outwall_temperature, fluxsolved, Sig, Re_function, Tcoolant, \
visccoolant, condcoolant, Cpmeth, rho, Vitesse, Pcoolant, LambdaTC, Celerite, \
hlnormal, singpertes, Pcoolant2 = mainsolver(
    Sig, rho, Tcoolant, visccoolant, condcoolant, Cpmeth,
    Pcoolant, LambdaTC, 1, reps_max + 1)

# Second and more iteration of the solving
# for i in range(0, reps_max):
#     newa = Sig[2]
#     Sig = [newa]
#     Tcoolant = []
#     visccoolant = []
#     condcoolant = []
#     Cpmeth = []
#     rho = []
#     newLambdatc = LambdaTC[1]
#     LambdaTC = []
#     Pcoolant = []
#     Tcoolant.append(Tl_init)
#     Pcoolant.append(Pl_init)
#     LambdaTC.append(newLambdatc)
#     visccoolant.append(meth.viscCH4(Pcoolant[0], Tcoolant[0], fluid))
#     condcoolant.append(meth.condCH4(Pcoolant[0], Tcoolant[0], fluid))
#     Cpmeth.append(meth.cpCH4(Pcoolant[0], Tcoolant[0], fluid))
#     rho.append(meth.rhoCH4(Pcoolant[0], Tcoolant[0], fluid))
#     hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, inwall_temperature, \
#     outwall_temperature, fluxsolved, Sig, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, \
#     Vitesse, Pcoolant, LambdaTC, Celerite, hlnormal, singpertes, Pcoolant2 = \
#         mainsolver(Sig, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, Pcoolant, LambdaTC, i + 2,
#                    reps_max + 1)

end_m = time.perf_counter()  # End of the main solution timer
time_elapsed_m = time.ctime(end_m - start_main_time)[14:19]  # Main elapsed time converted in minutes:secondes
start_d2 = time.perf_counter()  # Start of the display of 2D timer

print("█                                                                          █")
print("█ Execution time for the initialisation :", time_elapsed_i, "                           █")
print("█                                                                          █")
print("█ Execution time for the computation of the main solution :", time_elapsed_m, "         █")
print("█                                                                          █")
# %% Display of the first results
"Display of the results"
print("█ Display of results in 2D                                                 █")

colormap = plt.cm.magma
inv = 0, 0, 0  # 1 means should be reversed
view3d(inv, xcanaux, ycanaux, inwall_temperature, colormap, "Wall temperature on the gas side (in K)", size2,
       limitation)

Cel03 = [x * 0.3 for x in Celerite]
"""
plt.figure(dpi=300)
plt.plot(xcanauxre, Re_function, color='blue')
plt.title("Reynolds number of the coolant as a function of the engine axis")
plt.show()
"""
plt.figure(dpi=300)
plt.plot(xcanaux, hlcor, color='blue', label='Hl corrected')
plt.plot(xcanaux, hlnormal, color='cyan', label='Hl')
plt.title("Convection coeff as a function of the engine axis")
plt.legend()
plt.show()
"""
plt.figure(dpi=300)
plt.plot(xcanauxre, visc_function, color='orangered')
plt.title("Gas viscosity as a function of the engine axis")
plt.show()
plt.figure(dpi=300)
plt.plot(xcanauxre, cp_function, color='orangered')
plt.title("Gas Cp as a function of the engine axis")
plt.show()
plt.figure(dpi=300)
plt.plot(xcanauxre, lamb_function, color='orangered')
plt.title("Gas conductivity (lambda) as a function of engine axis")
plt.show()
plt.figure(dpi=300)
plt.plot(xcanauxre, Prandtl_function, color='orangered')
plt.title("Gas Prandtl number as a function of engine axis")
plt.show()
plt.figure(dpi=300)
plt.plot(xcanauxre, hg_function, color='orangered')
plt.title('Convection coefficient Hg as a function of engine axis')
plt.show()
plt.figure(dpi=300)
Sig.pop()
plt.plot(xcanauxre, Sig, color='orangered')
plt.title("Sigma as a function of the engine axis")
plt.show()
"""

plt.figure(dpi=300)
plt.plot(xcanaux, inwall_temperature[1:], color='magenta', label='Twg')
plt.plot(xcanaux, outwall_temperature[1:], color='orangered', label='Twl')
plt.title('Wall temperature (in K) as a function of engine axis')
plt.legend()
plt.show()

plt.figure(dpi=300)
Tcoolant.pop()
plt.plot(xcanaux, Tcoolant, color='blue')
plt.title('Coolant temperature (in K) as a function of engine axis')
plt.show()

"""
plt.figure(dpi=300)
rho.pop()
plt.plot(xcanauxre, rho, color='blue')
plt.title('Density of the coolant as a function of engine axis')
plt.show()
plt.figure(dpi=300)
visccoolant.pop()
plt.plot(xcanauxre, visccoolant, color='blue')
plt.title('Viscosity of the coolant as a function of engine axis')
plt.show()
plt.figure(dpi=300)
condcoolant.pop()
plt.plot(xcanauxre, condcoolant, color='blue')
plt.title('Conductivity of the coolant as a function of engine axis')
plt.show()
plt.figure(dpi=300)
Cpmeth.pop()
plt.plot(xcanauxre, Cpmeth, color='blue')
plt.title('Cp of the coolant as a function of engine axis')
plt.show()
"""

plt.figure(dpi=300)
plt.plot(xcanaux, Vitesse, color='blue', label='Coolant')
plt.plot(xcanaux, Cel03, color='orange', label='Mach 0.3 limit')
plt.title('Velocity of the coolant as a function of engine axis')
plt.legend()
plt.show()
Pcoolant.pop()
Pcoolant2.pop()
singpertes.pop()
plt.figure(dpi=300)
plt.plot(xcanaux, Pcoolant, color='orange', label='Pressure losses')
plt.title('Coolant pressure as a function of engine axis')
plt.legend()
plt.show()
LambdaTC.pop()
plt.figure(dpi=300)
plt.plot(xcanaux, LambdaTC, color='orangered')
plt.title('Thermal conductivity of the wall as a function of engine axis')
plt.show()
plt.figure(dpi=300)
plt.plot(xcanaux, Celerite, color='pink')
plt.title('Sound velocity of the coolant as a function of engine axis')
plt.show()

colormap = plt.cm.plasma
inv = 0, 0, 0
view3d(inv, xcanaux, ycanaux, fluxsolved, colormap, "Heat flux (in MW/m²)", size2, limitation)
colormap = plt.cm.coolwarm
inv = 0, 0, 0
view3d(inv, xcanaux, ycanaux, Tcoolant, colormap, "Temperature of the coolant", size2, limitation)

# %% Flux computation in 2D and 3D
"Computation for 2D graphs"
# 3D and 2D temperature contour plot
# At the beginning of the chamber
larg_ailette.reverse()
pas = larg_ailette[-1] + larg_canal[-1]
epaisseur = e_conv
hauteur = ht_canal[-1]
largeur = larg_canal[-1]
Hg = hg_function[-1]
Tg = hotgas_temp_list[-1]
Hl = hlnormal[-1]
Tl = Tcoolant[-1]
dx = 0.00004  # *3.5
lamb = LambdaTC[-1]
print("█                                                                          █")
print("█ Results at the beginning of the chamber :                                █")
where = " at the beginning of the chamber"
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 5, 1, 1, where)
# At the throat
pos_col = ycanaux.index(min(ycanaux))
pas = larg_ailette[pos_col] + larg_canal[pos_col]
epaisseur = e_col
hauteur = ht_canal[pos_col]
largeur = larg_canal[pos_col]
Hg = hg_function[pos_col]
Tg = hotgas_temp_list[pos_col]
Hl = hlnormal[pos_col]
Tl = Tcoolant[pos_col]
dx = 0.000025  # *3.5
lamb = LambdaTC[pos_col]
print("█                                                                          █")
print("█ Results at the throat :                                                  █")
where = " at the throat"
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 15, 1, 2, where)
# At the end of the divergent
pas = larg_ailette[0] + larg_canal[0]
epaisseur = e_tore
hauteur = ht_canal[0]
largeur = larg_canal[0]
Hg = hg_function[0]
Tg = hotgas_temp_list[0]
Hl = hlnormal[0]
Tl = Tcoolant[0]
dx = 0.00004
lamb = LambdaTC[0]
print("█                                                                          █")
print("█ Results at the end of the divergent :                                    █")
where = " at the end of the divergent"
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 5, 1, 1, where)

end_d2 = time.perf_counter()  # End of the display of 2D timer
time_elapsed_d2 = time.ctime(end_d2 - start_d2)[14:19]  # Display of 2D elapsed time converted in minutes:secondes

print("█                                                                          █")
print("█ Execution time for the display of 2D :", time_elapsed_d2, "                            █")
print("█                                                                          █")

"Computation for 3D graph"
choix = int(input("█ 3D temperature contour visualisation ? (1 = yes, 2 = no)                 █"))

if choix == 1:
    # 3D display
    start_d3 = time.perf_counter()
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

            lamb = LambdaTC[i]
            t3d = carto2D(larg_ailette[i] + larg_canal[i], wall_thickness[i], ht_canal[i], larg_canal[i], dx,
                          hg_function[i], lamb,
                          hotgas_temp_list[i], hlnormal[i], Tcoolant[i], 3, 0, 1, "")
            eachT.append(t3d)
            progressbar.update(1)

    carto3d([0, 0, 0], xcanaux, ycanaux, eachT, plt.cm.Spectral_r, '3D view of wall temperatures', nbc, limitation)

    end_d3 = time.perf_counter()  # End of the display of 3D timer
    time_elapsed_d3 = time.ctime(end_d3 - start_d2)[14:19]  # Display of 3D elapsed time converted in minutes:secondes
    print("█                                                                          █")
    print("█ Execution time for the display of 3D :", time_elapsed_d3, "                            █")
start_e = time.perf_counter()  # Start of the end timer
# %% Reversion of the different lists
print("█                                                                          █")

cross_section_area_list.reverse()
gamma_list.reverse()
mach_list.reverse()
hotgas_temp_list.reverse()
xcanaux.reverse()
ycanaux.reverse()
larg_canal.reverse()
ht_canal.reverse()
area_channel.reverse()
visc_function.reverse()
cp_function.reverse()
lamb_function.reverse()
hg_function.reverse()
Prandtl_function.reverse()
Sig.reverse()
inwall_temperature.reverse()
outwall_temperature.reverse()
fluxsolved.reverse()
Tcoolant.reverse()
Vitesse.reverse()
Re_function.reverse()
hlnormal.reverse()
rho.reverse()
visccoolant.reverse()
condcoolant.reverse()
Cpmeth.reverse()
Pcoolant.reverse()

# %% Preparation of the lists for CAD modelisation
"Changing the coordinates of the height of the channels (otherwise it is geometrically wrong)"
print("█                                                                          █")
angles = []
newxhtre = []
newyhtre = []

with tqdm(total=nb_points_channel,
          desc="█ Computing channel height     ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
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
        progressbar.update(1)

# Checking the height of channels
verification = []
with tqdm(total=nb_points_channel,
          desc="█ Checking channel height      ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points_channel):
        verifhtre = (((newxhtre[i] - xcanaux[i]) ** 2) + ((newyhtre[i] - ycanaux[i]) ** 2)) ** 0.5
        verification.append(verifhtre)
        progressbar.update(1)

plt.figure(dpi=300)
plt.plot(newxhtre, newyhtre, color='blue', label='New height')
plt.plot(xcanaux, ycanaux, color='chocolate', label='Former height')
plt.title("Geometrical aspect of the channel (height as a function of the engine axis)")
plt.axis("equal")
plt.legend()
plt.show()
plt.figure(dpi=300)
plt.plot(xcanaux, verification)
plt.title("Checking the height of the generated channels")
plt.show()

# %% Writing the results of the study in a CSV file

valuexport = open("valuexport.csv", "w", newline="")
geometry1 = open("geometry1.csv", "w", newline="")
geometry2 = open("geometry2.csv", "w", newline="")
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

with tqdm(total=nb_points,
          desc="█ Writing results in .csv      ",
          unit="|   █", bar_format="{l_bar}{bar}{unit}",
          ncols=76) as progressbar:
    for i in range(0, nb_points):
        if i < nb_points_channel:
            geometry1_writer.writerow((newxhtre[i] * (-1000), newyhtre[i] * 1000))
            geometry2_writer.writerow((ycanaux[i] * 1000, newxhtre[i] * (-1000)))
            valuexport_writer.writerow(
                (x_coord_list[i], y_coord_list[i], aire_saved[i], gamma_saved[i],
                 mach_list_saved[i], pressure_list[i],
                 hotgas_temperature_saved[i], xcanaux[i], ycanaux[i],
                 larg_canal[i], ht_canal[i], area_channel[i], visc_function[i],
                 cp_function[i], lamb_function[i], Prandtl_function[i],
                 hg_function[i], Sig[i], inwall_temperature[i],
                 outwall_temperature[i], fluxsolved[i], Tcoolant[i],
                 Re_function[i], hlnormal[i], rho[i], visccoolant[i],
                 condcoolant[i], Cpmeth[i], Vitesse[i], Pcoolant[i],
                 LambdaTC[i], newxhtre[i], newyhtre[i]))
        else:
            valuexport_writer.writerow(
                (x_coord_list[i], y_coord_list[i], aire_saved[i], gamma_saved[i],
                 mach_list_saved[i], pressure_list[i],
                 hotgas_temperature_saved[i], ' ', ' ', ' ', ' ', ' ', ' ',
                 ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                 ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))
        progressbar.update(1)

valuexport.close()
geometry1.close()
geometry2.close()

end_t = time.perf_counter()  # End of the total timer
time_elapsed_e = time.ctime(end_t - start_e)[14:19]  # End elapsed time converted in minutes:secondes
print("█                                                                          █")
print("█ Execution time of the end of the program :", time_elapsed_e, "                        █")

if choix == 1:
    time_elapsed_t_w3D = time.ctime((end_t - start_d3) + (end_d2 - start_time))[14:19]
    # Total elapsed time with 3D computation (without time waited to choose 3D) converted in minutes:secondes
    print("█                                                                          █")
    print("█ Total execution time with 3D computation :", time_elapsed_t_w3D, "                        █")

time_elapsed_t = time.ctime((end_t - start_e) + (end_d2 - start_time))[14:19]
# Total elapsed time (without time waited to choose 3D) converted in minutes:secondes

print("█                                                                          █")
print("█ Total execution time without 3D computation :", time_elapsed_t, "                     █")
print("█                                                                          █")
print("███████████████████████████████████ END ████████████████████████████████████")

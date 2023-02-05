# %% CTE
"""
Created on Fri Nov 27 14:47:27 2020

Original author: Julien S

Refactored and improved by Mehdi D, Paul B, Paul M, Eve X and Antoine R

WARNING: This Python file was rewritten only for the Viserion_2023 project.
Any changes might affect the results.
"""

# %% Import of packages
# # Calculations
# from math import *
# from sympy import Symbol, nsolve, Eq
# import sympy as mp
# import numpy as np
# from decimal import *  # ?
#
# # Graphics
# import matplotlib.pyplot as plt
# from pylab import *
# from mpl_toolkits.mplot3d import Axes3D
#
# # Files in the CTE folder
# # Solving
# from machsolve import Mach_solv
# from pressuresolve import Pressure_solv
# from temperaturesolve import Temperature_solv
# from musolve import mu
# from Tcorsolve import tempcorrige
# from heatequationsolve import *
# # from PrandtlVonKarman import PVK
# # Data
# from methan import *
# from Canaux import canauxangl, canaux
# # Graphics
# from ProgressBar import *
# from volume3d import *
# from graphic3d import view3d
#
# # Chemical species data
# from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
# import CoolProp.CoolProp as CoolProp
#
# # Others
# import csv  # To read .txt

import csv
import matplotlib.pyplot as plt
from pylab import *
from math import *
from machsolve import Mach_solv
from pressuresolve import Pressure_solv
from temperaturesolve import Temperature_solv
from musolve import mu
from Tcorsolve import tempcorrige
from ProgressBar import *
from Canaux import canaux
from sympy import Symbol, nsolve
import sympy as mp
from methan import *
from graphic3d import view3d
from heatequationsolve import *
from volume3d import *
from CoolProp.CoolProp import PropsSI

print("██████████████████████████ Cool The Engine V 2.0.0 █████████████████████████")
print("█                  Innovative Propulsion Laboratory - IPL                  █")
print("█                                                                          █")
print("█                                                                          █")
# %% Engine initialisation
Bar = ProgressBar(100, 30, "Initialization                  ")

"Viserion settings"

mesh_size = 0.25  # Distance between two points of calculation

x_coords_filename = f"input/{mesh_size}/x.txt"  # X coordinates of the Viserion
y_coords_filename = f"input/{mesh_size}/y.txt"  # Y coordinates of the Viserion
input_CEA_data = "input/Viserion_2023.txt"  # Viserion's parameters (found with CEA)

"Constant input_data_list"
lim22 = 600  # Not used anymore
size2 = 18  # Used for the height of the display in 3D view
ange = -18.884  # Not used anymore ?
lim11 = 200  # Not used anymore ?
epso = 10  # Not used anymore ?
machtype = 0  # 0 for one equation and else for another equation in calculation of mach
limitation = 0.05  # used to build the scales in 3D view

# %% Reading Viserion_2023.txt
crinit = csv.reader(open(input_CEA_data, "r"))
input_data_list = [row[1] for row in crinit]

# Store CEA output in lists
c_init = float(input_data_list[0])  # Sound velocity in the chamber
c_col = float(input_data_list[1])  # Sound velocity in the throat
debit_LOX = float(input_data_list[2])  # LOX debit
debit_LCH4 = float(input_data_list[3])  # CH4 debit
rho_init = float(input_data_list[4])  # Initial density of the gases
Pc = float(input_data_list[5])  # Pressure in the chamber
Tc = float(input_data_list[6])  # Combustion temperature (in the chamber?)
gamma_c = float(input_data_list[7])  # Gamma in the chamber
gamma_t = float(input_data_list[8])  # Gamma in the throat
gamma_e = float(input_data_list[9])  # Gamma at the exit
M = float(input_data_list[10])  # Molar mass of the gases
Cstar = float(input_data_list[11])  # Caracteristic velocity

# Store input dimensions in lists
Dcol = float(input_data_list[12])  # Convergent radius of curvature
Rcol = float(input_data_list[13])  # Throat radius of curvature
Ac = float(input_data_list[14])  # Throat diameter
DiamCol = float(input_data_list[15])  # Throat diameter

Bar.update(100)
print()
print("█                                                                          █")

# %% Import of the (X,Y) coordinates of the Viserion
Bar = ProgressBar(100, 30, "Import of (X,Y) coordinates     ")

# Reading the coordinate files
crx = csv.reader(open(x_coords_filename, "r"))
cry = csv.reader(open(y_coords_filename, "r"))

# Storing the X,Y coordinates in lists
x_value = [float(row[0]) / 1000 for row in crx]
y_value = [float(row[0]) / 1000 for row in cry]

Bar.update(100)
print()
print("█                                                                          █")

# Plot of the upper profile of the engine
"""
plt.figure(dpi=200)
plt.plot(x_value, y_value, color='black')
plt.title('Profile of the Viserion', color='black')
plt.show()
"""

# Computation and plot of the mesh density of the engine
"""
R = [abs(x_value[i] - x_value[i + 1]) for i in range(0, len(x_value) - 1)]
R.append(R[-1])
colooo = plt.cm.binary
inv = 1, 1, 1  # 1 means should be reversed
view3d(inv, x_value, y_value, R, colooo, 'Mesh density', size2 - 2, limitation)
"""

# %% Computation of the cross-sectional areas of the engine
Bar = ProgressBar(100, 30, "Sectionnal areas computation    ")

aire = [pi * r ** 2 for r in y_value]

Bar.update(100)
print()

# Plots of the cross-sectionnal areas
"""
plt.figure(dpi=200)
plt.plot(x_value,aire,color='black')
plt.title("Area in m² as a function of engine axis")
plt.show()
"""
print("█                                                                          █")

# %% Adiabatic constant (gamma) parametrization
"Linear interpolation of gamma"
Bar = ProgressBar(100, 30, "Gamma computation               ")

i = 0  # Index of the beginning of the convergent
a = 1
b = 1
while a == b:  # Read y values two per two in order to detect the beginning of the convergent
    a = y_value[i]
    i += 1
    b = y_value[i]
# Gamma in the cylindrical chamber
gamma = []
for i_throat in range(0, i):  # Gamma is constant before the beginning of the convergent along the chamber
    gamma.append(gamma_c)

# Gamma in the convergent
i_throat = y_value.index(min(y_value))  # Throat index
k = i_throat - i  # Number of points in the convergent
c = gamma_c
for m in range(-1, k - 1):
    # Linear interpolation between beginning and end of convergent:
    # (yi+1)=((y2-y1)/(x2-x1))*abs((xi+1)-(xi))
    l = ((gamma_t - gamma_c) / (x_value[i_throat] - x_value[i])) * abs(x_value[i + 1 + m] - x_value[i + m])
    c += l
    gamma.append(c)

# Gamma in the divergent nozzle
p = len(x_value) - i_throat  # Number of points in the divergent
t = gamma_t
for q in range(-1, p - 1):  # Linear interpolation between beginning and end of divergent
    n = ((gamma_e - gamma_t) / (x_value[-1] - x_value[i_throat])) * abs(
        x_value[i_throat + 1 + q] - x_value[i_throat + q])
    t += n
    gamma.append(t)

Bar.update(100)
print()

# Plot of the gamma linearisation
"""
#print(gamma)
#print(len(gamma))
plt.figure(dpi=200)
plt.plot(x_value,gamma,color='gold')
plt.title("Gamma linearisation")
plt.show()
"""
print("█                                                                          █")

# %% Mach number computation
"Computation of the initial velocity and mach number of the gases"
v_init = (debit_LOX + debit_LCH4) / (rho_init * aire[0])  # Initial velocity of the gases
M_init = v_init / c_init  # Initial mach number
M1 = M_init
mach_function = [M_init]
b = 0
Bar = ProgressBar(100, 30, "Mach number computation         ")
long = len(x_value)  # Number of points
av = 100 / (long - 1)

"Mach number computations along the engine"
for i in range(0, long - 1):
    A1 = aire[i]
    A2 = aire[i + 1]
    pos = i
    g = Mach_solv(A1, A2, M1, gamma[i], pos, machtype)
    mach_function.append(g)
    M1 = g
    b += av
    Bar.update(b)

print()

# Plots of the Mach number in the engine (2D/3D)
plt.figure(dpi=200)
plt.plot(x_value, mach_function, color='gold')
plt.title("Mach number as a function of the engine axis")
plt.show()

colooo = plt.cm.Spectral
inv = 1, 1, 1  # 1 means should be reversed
view3d(inv, x_value, y_value, mach_function, colooo, 'Mach number', size2 - 2, limitation)

print("█                                                                          █")

# %% Static pressure computation
"Static pressure computation"
c = 0
Bar = ProgressBar(100, 30, "Static pressure computation     ")
ac = 100 / (long - 1)

pressure_function = [Pc]  # (in Pa)
"Static pressure computations along the engine"
for i in range(0, long - 1):
    if i == long + 1:
        M1 = mach_function[i]
        M2 = mach_function[i]
    else:
        M1 = mach_function[i]
        M2 = mach_function[i + 1]
    P1 = pressure_function[i]
    P2 = Pressure_solv(M1, M2, P1, gamma[i])
    pressure_function.append(P2)
    c += ac
    Bar.update(c)

print()

# Plot of the static pressure (2D/3D)
plt.figure(dpi=200)
plt.plot(x_value, pressure_function, color='gold')
plt.title("Static pressure as a function of the engine axis")
plt.show()

colooo = plt.cm.gist_rainbow_r
inv = 1, 1, 1  # 1 means should be reversed
view3d(inv, x_value, y_value, pressure_function, colooo, 'Static pressure', size2 - 2, limitation)

print("█                                                                          █")

# %% Temperature computation
b = 0
Bar = ProgressBar(100, 30, "Temperature computation         ")
ay = 100 / (long - 1)

# Hot gas temperature computations along the engine
hotgas_temperature = [Tc]
for i in range(0, long - 1):
    if i == long + 1:
        M1 = mach_function[i]
        M2 = mach_function[i]
    else:
        M1 = mach_function[i]
        M2 = mach_function[i + 1]
    T1 = hotgas_temperature[i]
    T2 = Temperature_solv(M1, M2, T1, gamma[i])
    hotgas_temperature.append(T2)
    b += ay
    Bar.update(b)

print()
# List of corrected gas temperatures (max diff with original is about 75 K)
hotgas_temp_corrected = [tempcorrige(hotgas_temperature[i], gamma[i], mach_function[i]) for i in range(0, long)]

# Plots of the temperature in the engine (2D/3D)
"""
plt.figure(dpi=200)
plt.plot(x_value,hotgas_temperature,color='gold')
plt.title("Temperature as a function of the engine axis")
plt.show()

colooo = plt.cm.terrain_r
inv = 1, 1, 1  # 1 means should be reversed
view3d(inv,x_value,y_value,hotgas_temperature,colooo,'Temperature of the gases',size2-2,limitation)
"""

print("█                                                                          █")

# %% Canal parameters
"Number of channels and tore position"
nbc = 40  # Number of channels
tore = 0.103  # Position of the manifol from the throat (in m)

"Width of the channels"
lrg_c2 = 0.003  # Width of the canal in at the injection plate (in m)
lrg_c = 0.003  # Width of the canal at the end of the cylindrical chamber (in m)
lrg_col = 0.002  # Width of the canal in the throat (in m)
lrg_div = 0.003  # Width of the canal at the extremity of the nozzle (in m)

"Height of the channels"
ht_c2 = 0.003  # Height of the canal at the injection plate (in m)
ht_c = 0.003  # Height of the canal at the end of the cylindrical chamber (in m)
ht = 0.002  # Height of the canal in the throat (in m)
ht_div = 0.003  # Height of the canal at the extremity of the nozzle (in m)

# %% Thicknesses
e_c = 0.001  # Thickness of the chamber (in m)
e_col = 0.001  # Thickness of the throat (in m)
e_div = 0.001  # Thickness of the divergent (in m)

# %% Growth factors
n1 = 1  # Width convergent
n2 = 1  # Width divergent
n3 = 1  # Height convergent
n4 = 1  # Height divergent
n5 = 1  # Thickness convergent
n6 = 1  # Thickness divergent

# %% Material
"Thermal conductivity of the material"
material = 0
if material == 0:  # Pure copper
    condcoeff1 = -0.065665283166
    condcoeff2 = 421.82710859
elif material == 1:  # CuCrZr
    condcoeff1 = -0.0269
    condcoeff2 = 365.33

# %% Coolant
"Properties of the coolant"
fluid = "Methane"
rho_initCH4 = 425  # Density of the CH4
If_reg = debit_LCH4  # Total mass flow rate (in kg/s)
Tl_init = 110  # Initial temperature of the coolant (in K)
debit_total = If_reg / rho_initCH4  # Total volumic flow rate of the coolant (in m3/s)
Pl_init = 3700000  # Initial pressure of the coolant (in Pa)
Ro = 3  # Roughness (micrometers)

# %% Computation of canals
Bar = ProgressBar(100, 30, "Canal geometric computation     ")

"""Methode 2"""
xcanauxre, ycanauxre, larg_canalre, Areare, htre, reste, epaiss_chemise \
    = canaux(x_coords_filename, y_coords_filename, nbc, lrg_col, lrg_c,
             lrg_div, ht, ht_c, ht_div, tore, debit_total, n1, n2, n3, n4,
             e_col, e_div, e_c, n5, n6, lrg_c2, ht_c2)

"""Methode 1"""
"""
xcanauxre, ycanauxre, larg_canalre, Areare, htre = canauxangl(x_coords_filename, y_coords_filename,
                                                              nbc, lrg_col, ht, ht_c, ht_div, tore,
                                                              debit_total, epaisseur_chemise)
"""

Bar.update(100)
print()
print("█                                                                          █")

# We reverse the data in order to calculate from the manifold to the injection (x is in reverse)
epaiss_chemise.reverse()
xcanauxre.reverse()
larg_canalre.reverse()
Areare.reverse()
htre.reverse()
ycanauxre.reverse()
fin = len(xcanauxre)  # Index of the end of the channels, after removing everything before the manifold

# Save the data for exporting, before altering the original lists
hotgas_temperature_saved = hotgas_temperature[:]
aire_saved = aire[:]
mach_function_saved = mach_function[:]
gamma_saved = gamma[:]

# Remove the data points before the manifold
hotgas_temperature_ = hotgas_temperature[:fin]
aire_ = aire[:fin]
mach_function_ = mach_function[:fin]
gamma_ = gamma[:fin]

gamma.reverse()
mach_function.reverse()
aire.reverse()
hotgas_temperature.reverse()


def mainsolver(Sig, b, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, ay, Pcoolant, LambdaTC, entropy):
    """
    This is the main function used for solving the 1D case.
    The geometry is discretised into a 1 dimensionnal set of points.
    The function uses a marching algorithm, and computes all the relevant physical
    quantities at each point. The values obtained are then used on the next point.
    """

    # Lists containing the physical quantities at each point
    hlcor = []
    visc_function = []
    cp_function = []
    Re_function = []
    lamb_function = []
    Prandtl_function = []
    hg_function = []
    inwall_temperature = []
    outwall_temperature = []
    fluxsolved = []
    Vitesse = []
    Celerite = []
    hlnormal = []
    error_D_ = []
    singpertes = [Pcoolant[0]]
    Pcoolant2 = [Pcoolant[0]]
    phase = 0
    positioncol = ycanauxre.index(min(ycanauxre))

    # Main computation loop
    for i in range(0, len(xcanauxre)):
        Lambda_tc = LambdaTC[i]
        x = xcanauxre[i]
        c = larg_canalre[i]

        # Hydraulic diameter (4*Area/Perimeter)
        Dhy = (2 * htre[i] * c) / (htre[i] + c)

        # Velocity of the coolant
        V = debit_LCH4 / (nbc * rho[i] * Areare[i])
        Vitesse.append(V)

        # Reynolds number
        Re = (V * Dhy * rho[i]) / visccoolant[i]
        Re_function.append(Re)

        # Prandtl number
        Pr_cool = (visccoolant[i] * Cpmeth[i]) / condcoolant[i]

        T1 = hotgas_temperature[i]
        # Compute viscosity, Cp, conductivity and Prandtl number with musolve.py
        viscosite, cp, lamb, Pr = mu(T1, M, gamma[i])

        # Store the results
        visc_function.append(viscosite)
        cp_function.append(cp)
        lamb_function.append(lamb)
        Prandtl_function.append(Pr)

        # Gas-side convective heat transfer coefficient (Bartz equation)
        hg = (0.026 / (DiamCol ** 0.2) * (((viscosite ** 0.2) * cp) / (Pr ** 0.6)) * (
                (Pc / Cstar) ** 0.8) * ((DiamCol / Dcol) ** 0.1) * ((Ac / aire[i]) ** 0.9)) * Sig[i]
        hg_function.append(hg)

        # Radiative heat flux assuming black body radiation
        Tg = hotgas_temperature[i]
        steff = 5.6697 * 10 ** (-8)  # Stefan-Boltzmann constant
        emissivity = 0.02  # 2% emissivity for CH4
        qr = emissivity * steff * (T1 ** 4)  # Radiative heat flux

        # # This is not used, didn't delete because it might have a purpose later
        # if i + 1 == len(xcanauxre):
        #     EtalArea = (2 * c + 2 * htre[i]) * abs(xcanauxre[i] - xcanauxre[i - 1])
        # else:
        #     EtalArea = (2 * c + 2 * htre[i]) * abs(xcanauxre[i + 1] - x)

        Gdeb = debit_LCH4 / (Areare[i] * nbc)

        # In case of single-phase flow
        if phase == 0 or phase == 2:
            f = (1.82 * log10(Re) - 1.64) ** (-2)
            Nu = ((f / 8) * (Re - 1000) * Pr_cool) / (1 + 12.7 * ((f / 8) ** 0.5) * ((Pr_cool ** (2 / 3)) - 1))
            hl = Nu * (condcoolant[i] / Dhy)
            Xqual = PropsSI("Q", "P", Pcoolant[i], "T", Tcoolant[i], fluid)

        # In case of two-phase flow (boiling etc.)
        else:
            # Grab all the fluid data using CoolProp
            Base0 = PropsSI("H", "Q", 0, "P", Pcoolant[i], fluid)
            Base1 = PropsSI("H", "Q", 1, "P", Pcoolant[i], fluid)
            Xqual = (H2 - Base0) / (Base1 - Base0)
            rhogp = PropsSI("D", "Q", 1, "P", Pcoolant[i], fluid)
            rholp = PropsSI("D", "Q", 0, "P", Pcoolant[i], fluid)
            mugp = PropsSI("V", "Q", 1, "P", Pcoolant[i], fluid)
            mulp = PropsSI("V", "Q", 0, "P", Pcoolant[i], fluid)
            Xeq = (((1 - Xqual) / Xqual) ** 0.9) * ((rhogp / rholp) ** 0.5) * ((mugp / mulp) ** 0.1)
            SurfaceT = PropsSI("I", "Q", 0, "P", Pcoolant[i], fluid)
            Kp = float(Pcoolant[i]) / ((SurfaceT * 9.81 * (rholp - rhogp)) ** 0.5)
            We = (Gdeb ** 2) * Dhy / (rholp * SurfaceT)
            # 1-0 H
            HVP = PropsSI("H", "Q", 1, "P", Pcoolant[i], fluid)
            HLP = PropsSI("H", "Q", 0, "P", Pcoolant[i], fluid)
            HHV = HVP - HLP

            flux1 = flux + 0.5
            flux2 = flux
            while (abs(flux1 - flux2) / flux2) > 0.000001:
                flux1 = (flux1 + flux2) / 2
                Bo = 1000000 * flux1 / (Gdeb * HHV)

                # Use a different correlation for the Nusselt number,
                # depending on the vapor quality
                if Xqual < 0.6:
                    Nu = 12.46 * (Bo ** 0.544) * (We ** 0.035) * (Kp ** 0.614) * (Xeq ** 0.031)
                else:
                    Nu = 0.00136 * (Bo ** (-1.442)) * (We ** 0.074)

                # Compute the convective heat-transfer coefficient
                hl = Nu * (condcoolant[i] / Dhy)

                # Compute dimensions of the fins
                D = 2 * (ycanauxre[i] - epaiss_chemise[i])
                d = (pi * (D + htre[i] + epaiss_chemise[i]) - nbc * c) / nbc
                m = ((2 * hl) / (d * Lambda_tc)) ** 0.5

                # Corrected coefficient, taking the fin effect into account
                hl_cor = hl * ((nbc * c) / (pi * D)) + nbc * (
                        (2 * hl * Lambda_tc * (((pi * D) / nbc) - c)) ** 0.5) * ((tanh(m * htre[i])) / (pi * D))

                # Save the data in lists
                hg = hg_function[i]
                hl = hl_cor
                Tl = Tcoolant[i]
                e = epaiss_chemise[i]
                L = Lambda_tc
                mp.dps = 150

                # Use sympy to solve a system of 2 equations and 2 unknowns
                cx1 = Symbol('cx1')
                cx2 = Symbol('cx2')
                f1 = hg * (Tg - cx1) - (L / e) * (cx1 - cx2)
                f2 = hl * (cx2 - Tl) - (L / e) * (cx1 - cx2)

                # Solve the system numerically, giving an initial guess
                x_, y_ = nsolve((f1, f2), (cx1, cx2), (900, 700))
                # Flow computation
                flux2 = hl * (y_ - Tcoolant[i]) * 0.000001

        # %% Computations

        # Compute coolant-side convective heat-transfer coefficient
        hl = Nu * (condcoolant[i] / Dhy)
        hlnormal.append(hl)

        # Fin dimensions
        D = 2 * (ycanauxre[i] - epaiss_chemise[i])
        d = (pi * (D + htre[i] + epaiss_chemise[i]) - nbc * c) / nbc
        m = ((2 * hl) / (d * Lambda_tc)) ** 0.5

        # Correct for the fin effect
        hl_cor = hl * ((nbc * c) / (pi * D)) + nbc * \
                 ((2 * hl * Lambda_tc * (((pi * D) / nbc) - c)) ** 0.5) * (
                         (tanh(m * htre[i])) / (pi * D))
        hlcor.append(hl_cor)

        # Store the results
        hg = hg_function[i]
        hl = hlcor[i]
        Tl = Tcoolant[i]
        e = epaiss_chemise[i]

        # Use sympy to solve a system of 2 equations to solve the coupled
        # heat transfer in the channels
        L = Lambda_tc
        mp.dps = 150
        cx1 = Symbol('cx1')
        cx2 = Symbol('cx2')
        f1 = hg * (Tg - cx1) - (L / e) * (cx1 - cx2)
        f2 = hl * (cx2 - Tl) - (L / e) * (cx1 - cx2)
        x_, y_ = nsolve((f1, f2), (cx1, cx2), (700, 500))

        # Temperature at the walls
        inwall_temperature.append(x_)
        outwall_temperature.append(y_)

        # Compute heat flux throught the coolant side
        flux = hl * (y_ - Tcoolant[i]) * 0.000001
        fluxsolved.append(flux)

        Tw = inwall_temperature[i]
        Ts = hotgas_temperature[positioncol]
        Mak = mach_function[i]
        sigm = 1 / ((((Tw / (2 * Ts)) * (1 + (((gamma[i] - 1) / 2) * (Mak ** 2))) + 0.5) ** 0.68) * (
                (1 + (((gamma[i] - 1) / 2) * (Mak ** 2))) ** 0.12))
        Sig.append(sigm)
        Lambdatc = (condcoeff1 * ((Tw + outwall_temperature[i]) * 0.5) + condcoeff2)
        LambdaTC.append(Lambdatc)
        if i == xcanauxre.index(xcanauxre[-1]):
            Distance = ((xcanauxre[i - 1] - xcanauxre[i]) ** 2 + (ycanauxre[i - 1] - ycanauxre[i]) ** 2) ** 0.5
            xa = Distance
            ya = (2 * pi * ycanauxre[i - 1]) / nbc
            za = (2 * pi * ycanauxre[i]) / nbc
            perim = (2 * pi * ((aire[i] / pi) ** 0.5)) / nbc
            dA = xa * (2 * c + 2 * htre[i])
            timer = (((xcanauxre[i - 1] - xcanauxre[i]) ** 2 + (ycanauxre[i - 1] - ycanauxre[i]) ** 2) ** 0.5) / V
        else:
            Distance = ((xcanauxre[i + 1] - xcanauxre[i]) ** 2 + (ycanauxre[i + 1] - ycanauxre[i]) ** 2) ** 0.5
            xa = Distance
            ya = (2 * pi * ycanauxre[i + 1]) / nbc
            za = (2 * pi * ycanauxre[i]) / nbc

            perim = (2 * pi * ((aire[i] / pi) ** 0.5)) / nbc
            # dA=xa*perim

            dA = xa * (2 * c + 2 * htre[i])
            timer = (((xcanauxre[i + 1] - xcanauxre[i]) ** 2 + (ycanauxre[i + 1] - ycanauxre[i]) ** 2) ** 0.5) / V
        Q = abs(flux) * 1000000 * abs(dA)
        density = rho[i]
        debitmass = debit_LCH4 / nbc
        Tfu = (Q / (debitmass * Cpmeth[i])) + Tcoolant[i]
        if phase == 0 or phase == 2:
            rho2 = PropsSI("D", "P", Pcoolant[i], "T", Tfu, fluid)
            if (rho[i] / rho2) > 1.5:
                rho2 = PropsSI("D", "Q", 0, "T", Tfu, fluid)
            Re_sp = Dhy * Gdeb / visccoolant[i]
            fsp1 = 1
            fsp2 = (1 / (-2 * log10(((Ro / Dhy) / 3.7) + 2.51 / (Re_sp * (fsp1 ** 0.5))))) ** 2
            while abs((fsp1 / fsp2) - 1) > 0.0000001:
                fsp1 = fsp2
                fsp2 = (1 / (-2 * log10(((Ro / Dhy) / 3.7) + 2.51 / (Re_sp * (fsp1 ** 0.5))))) ** 2
            DeltaPfr = fsp1 * (Gdeb ** 2) * Distance / (2 * rho[i] * Dhy)
            DeltaPac = ((Gdeb ** 2) / (2 * rho[i])) * ((rho[i] / rho2) - 1)
        else:
            Re_tp = Dhy * Gdeb / visccoolant[i]
            eps_prim = (rholp / rhogp) / ((1 / Xqual) + ((rholp / rhogp) - 1))
            DeltaPac = (Gdeb ** 2) * Distance * (
                    (((1 - Xqual) ** 2) / (rholp * (1 - eps_prim))) + ((Xqual ** 2) / (rhogp * eps_prim)))
            rhotp = (rhogp * rholp) / (rhogp * (1 - Xqual) + rholp * Xqual)
            fsp1 = 1
            fsp2 = (1 / (-2 * log10(((Ro / Dhy) / 3.7) + (2.51 / (Re_tp * (fsp1 ** 0.5)))))) ** 2
            while abs((fsp1 / fsp2) - 1) > 0.0000001:
                fsp1 = fsp2
                fsp2 = (1 / (-2 * log10(((Ro / Dhy) / 3.7) + (2.51 / (Re_tp * (fsp1 ** 0.5)))))) ** 2
            DeltaPfr = fsp1 * (Gdeb ** 2) * Distance / (2 * rhotp * Dhy)
        NewPressure = Pcoolant[i] - (DeltaPfr + DeltaPac)
        Pcoolant.append(NewPressure)
        # Phase verification of the coolant at the next point using the raise in enthalpy
        if phase == 0:
            H1 = PropsSI("H", "P", Pcoolant[i], "T", Tcoolant[i], fluid)
            H2 = H1 + Q / debitmass
            Quality = PropsSI("Q", "H", H2, "P", Pcoolant[i + 1], fluid)
            if 0 < Quality < 1:
                phase = 1
            else:
                phase = 0
        else:
            H1 = H2
            H2 = H1 + Q / debitmass
            Quality = PropsSI("Q", "H", H2, "P", Pcoolant[i + 1], fluid)
            if 0 < Quality < 1:
                phase = 1
            else:
                phase = 2
        if phase == 0 or phase == 2:
            Tcoolant.append(Tfu)
        else:
            Tcoolant.append(Tcoolant[i])
        # CH4 new properties computation
        if phase == 0 or phase == 2:
            newvisc = PropsSI("V", "P", Pcoolant[i + 1], "T", Tcoolant[i], fluid)
            visccoolant.append(newvisc)
            newcond = PropsSI("L", "P", Pcoolant[i + 1], "T", Tcoolant[i], fluid)
            condcoolant.append(newcond)
            newcp = PropsSI("C", "P", Pcoolant[i + 1], "T", Tcoolant[i], fluid)
            Cpmeth.append(newcp)
            dens = PropsSI("D", "P", Pcoolant[i + 1], "T", Tcoolant[i], fluid)
            rho.append(dens)
            cele = PropsSI("A", "P", Pcoolant[i + 1], "T", Tcoolant[i], fluid)
            Celerite.append(cele)
        else:
            newvisc1 = PropsSI("V", "Q", 1, "P", Pcoolant[i + 1], fluid)
            newcond1 = PropsSI("L", "Q", 1, "P", Pcoolant[i + 1], fluid)
            newcp1 = PropsSI("C", "Q", 1, "P", Pcoolant[i + 1], fluid)
            dens1 = PropsSI("D", "Q", 1, "P", Pcoolant[i + 1], fluid)
            newvisc2 = PropsSI("V", "Q", 0, "P", Pcoolant[i + 1], fluid)
            newcond2 = PropsSI("L", "Q", 0, "P", Pcoolant[i + 1], fluid)
            newcp2 = PropsSI("C", "Q", 0, "P", Pcoolant[i + 1], fluid)
            dens2 = PropsSI("D", "Q", 0, "P", Pcoolant[i + 1], fluid)
            visccoolant.append(Quality * newvisc1 + (1 - Quality) * newvisc2)
            condcoolant.append(Quality * newcond1 + (1 - Quality) * newcond2)
            Cpmeth.append(Quality * newcp1 + (1 - Quality) * newcp2)
            rho.append(Quality * dens1 + (1 - Quality) * dens2)
            cele1 = PropsSI("A", "Q", 1, "P", Pcoolant[i + 1], fluid)
            cele2 = PropsSI("A", "Q", 0, "P", Pcoolant[i + 1], fluid)
            cele = Quality * cele1 + (1 - Quality) * cele2
            Celerite.append(cele)

        b += ay
        Bar.update(b)

    return hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, inwall_temperature, \
           outwall_temperature, fluxsolved, Sig, b, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, \
           Vitesse, Pcoolant, LambdaTC, Celerite, hlnormal, error_D_, singpertes, Pcoolant2


Sig = []
Tcoolant = []
Pcoolant = []
visccoolant = []
condcoolant = []
Cpmeth = []
rho = []
b = 0
Sig.append(1)
Tcoolant.append(Tl_init)
Pcoolant.append(Pl_init)
visccoolant.append(viscCH4(Pcoolant[0], Tcoolant[0], fluid))
condcoolant.append(condCH4(Pcoolant[0], Tcoolant[0], fluid))
Cpmeth.append(cpCH4(Pcoolant[0], Tcoolant[0], fluid))
rho.append(rhoCH4(Pcoolant[0], Tcoolant[0], fluid))
LambdaTC = [330]
entropy = [entCH4(Pcoolant[0], Tcoolant[0], fluid)]
b = 0
rep = 2
Bar = ProgressBar(100, 30, "Résolution globale              ")
ay = 100 / ((1 + rep) * len(xcanauxre))
hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, inwall_temperature, \
outwall_temperature, fluxsolved, Sig, b, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, Vitesse, \
Pcoolant, LambdaTC, Celerite, hlnormal, error_D_, singpertes, Pcoolant2 = mainsolver(
    Sig, b, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, ay, Pcoolant, LambdaTC, entropy)
# Second itération of the solving
for i in range(0, rep):
    newa = Sig[2]
    Sig = []
    Tcoolant = []
    visccoolant = []
    condcoolant = []
    Cpmeth = []
    rho = []
    newLambdatc = LambdaTC[1]
    LambdaTC = []
    Sig.append(newa)
    Pcoolant = []
    Tcoolant.append(Tl_init)
    Pcoolant.append(Pl_init)
    LambdaTC.append(newLambdatc)
    visccoolant.append(viscCH4(Pcoolant[0], Tcoolant[0], fluid))
    condcoolant.append(condCH4(Pcoolant[0], Tcoolant[0], fluid))
    Cpmeth.append(cpCH4(Pcoolant[0], Tcoolant[0], fluid))
    rho.append(rhoCH4(Pcoolant[0], Tcoolant[0], fluid))
    entropy = [entCH4(Pcoolant[0], Tcoolant[0], fluid)]
    hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, inwall_temperature, \
    outwall_temperature, fluxsolved, Sig, b, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, \
    Vitesse, Pcoolant, LambdaTC, Celerite, hlnormal, error_D_, singpertes, Pcoolant2 = \
        mainsolver(Sig, b, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, ay, Pcoolant, LambdaTC, entropy)

print()
print("█                                                                          █")

# %% Display of the first results
"Display of the results"
print("█ Display of results in 2D :                                               █")
print("█                                                                          █")

colooo = plt.cm.magma
inv = 0, 0, 0  # 1 means should be reversed
view3d(inv, xcanauxre, ycanauxre, inwall_temperature, colooo, "Wall temperature on the gas side", size2, limitation)

Cel03 = [x * 0.3 for x in Celerite]
"""
plt.figure(dpi=200)
plt.plot(xcanauxre, Re_function, color='blue')
plt.title("Reynolds number as a function of the engine axis")
plt.show()
"""
plt.figure(dpi=200)
plt.plot(xcanauxre, hlcor, color='blue', label='Hl corrigé')
plt.plot(xcanauxre, hlnormal, color='cyan', label='Hl')
plt.title("Convection coefficient Hl as a function of the engine axis")
plt.legend()
plt.show()
"""
plt.figure(dpi=200)
plt.plot(xcanauxre, visc_function, color='orangered')
plt.title("Gas viscosity as a function of the engine axis")
plt.show()
plt.figure(dpi=200)
plt.plot(xcanauxre, cp_function, color='orangered')
plt.title("Gas Cp as a function of the engine axis")
plt.show()
plt.figure(dpi=200)
plt.plot(xcanauxre, lamb_function, color='orangered')
plt.title("Gas conductivity (lambda)")
plt.show()
plt.figure(dpi=200)
plt.plot(xcanauxre, Prandtl_function, color='orangered')
plt.title("Prandtl number of gas")
plt.show()
plt.figure(dpi=200)
plt.plot(xcanauxre, hg_function, color='orangered')
plt.title('Convection coefficient Hg')
plt.show()
plt.figure(dpi=200)
Sig.pop()
plt.plot(xcanauxre, Sig, color='orangered')
plt.title("Sigma as a function of the engine axis")
plt.show()
"""

plt.figure(dpi=200)
plt.plot(xcanauxre, inwall_temperature, color='magenta', label='Twg')
plt.plot(xcanauxre, outwall_temperature, color='orangered', label='Twl')
plt.title('Wall temperature')
plt.legend()
plt.show()

plt.figure(dpi=200)
Tcoolant.pop()
plt.plot(xcanauxre, Tcoolant, color='blue')
plt.title('Coolant temperature')
plt.show()

"""plt.figure(dpi=200)
rho.pop()
plt.plot(xcanauxre, rho, color='blue')
plt.title('Volumic mass of the coolant')
plt.show()
plt.figure(dpi=200)
visccoolant.pop()
plt.plot(xcanauxre, visccoolant, color='blue')
plt.title('Viscosity of the coolant')
plt.show()
plt.figure(dpi=200)
condcoolant.pop()
plt.plot(xcanauxre, condcoolant, color='blue')
plt.title('Conductivity of the coolant')
plt.show()
plt.figure(dpi=200)
Cpmeth.pop()
plt.plot(xcanauxre, Cpmeth, color='blue')
plt.title('Cp of the coolant')
plt.show()"""

plt.figure(dpi=200)
plt.plot(xcanauxre, Vitesse, color='blue', label='Coolant')
plt.plot(xcanauxre, Cel03, color='orange', label='Mach 0.3 limit')
plt.title('Velocity of the coolant')
plt.legend()
plt.show()
Pcoolant.pop()
Pcoolant2.pop()
singpertes.pop()
plt.figure(dpi=200)
plt.plot(xcanauxre, Pcoolant, color='orange', label='Pressure losses')
plt.title('Pressure of the coolant')
plt.legend()
plt.show()
LambdaTC.pop()
plt.figure(dpi=200)
plt.plot(xcanauxre, LambdaTC, color='orangered')
plt.title('Thermal conductivity of the wall')
plt.show()
plt.figure(dpi=200)
plt.plot(xcanauxre, Celerite, color='pink')
plt.title('Sound velocity of the coolant')
plt.show()

colooo = plt.cm.plasma
inv = 0, 0, 0
view3d(inv, xcanauxre, ycanauxre, fluxsolved, colooo, "Heat flux (in MW/m²)", size2, limitation)
colooo = plt.cm.coolwarm
inv = 0, 0, 0
view3d(inv, xcanauxre, ycanauxre, Tcoolant, colooo, "Temperature of the coolant", size2, limitation)

# %% Flux computation in 2D and 3D
"Computation for 2D graphs"
# 3D and 2D cartography establishement of temperatures
# At the beginning of the chamber
reste.reverse()
pas = reste[-1] + larg_canalre[-1]
epaisseur = e_c
hauteur = htre[-1]
largeur = larg_canalre[-1]
Hg = hg_function[-1]
Tg = hotgas_temperature[-1]
Hl = hlnormal[-1]
Tl = Tcoolant[-1]
dx = 0.00004  # *3.5
lamb = LambdaTC[-1]
print("█ Results at the beginning of the chamber :                                █")
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 5, 1, 1)
# At the throat
poscol = ycanauxre.index(min(ycanauxre))
pas = reste[poscol] + larg_canalre[poscol]
epaisseur = e_col
hauteur = htre[poscol]
largeur = larg_canalre[poscol]
Hg = hg_function[poscol]
Tg = hotgas_temperature[poscol]
Hl = hlnormal[poscol]
Tl = Tcoolant[poscol]
dx = 0.000025  # *3.5
lamb = LambdaTC[poscol]
print("█ Results at the throat :                                                  █")
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 15, 1, 2)
# At the end of the divergent
pas = reste[0] + larg_canalre[0]
epaisseur = e_div
hauteur = htre[0]
largeur = larg_canalre[0]
Hg = hg_function[0]
Tg = hotgas_temperature[0]
Hl = hlnormal[0]
Tl = Tcoolant[0]
dx = 0.00004
lamb = LambdaTC[0]
print("█ Results at the end of the divergent :                                    █")
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 5, 1, 1)

"Computation for 3D graph"
choix2 = int(input("█ Visualisation des cartographie 3D ? (1=oui sinon 2)                      █"))
if choix2 == 1:
    # 3D display
    eachT = []
    for i in range(0, len(xcanauxre)):
        # print(i)
        lim1 = 0
        lim2 = 650
        pas = reste[i] + larg_canalre[i]
        epaisseur = epaiss_chemise[i]
        hauteur = htre[i]
        largeur = larg_canalre[i]
        Hg = hg_function[i]
        Tg = hotgas_temperature[i]
        Hl = hlnormal[i]
        Tl = Tcoolant[i]
        if lim2 >= i >= lim1:
            dx = 0.0001
        else:
            dx = 0.0001
        lamb = LambdaTC[i]
        t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 3, 0, 1)
        eachT.append(t3d)
        inv = [0, 0, 0]
        x = xcanauxre
    yprim = ycanauxre
    temp = eachT
    colooo = plt.cm.Spectral_r
    title = '3D view of wall temperatures'
    number = nbc
    carto3d(inv, x, yprim, temp, colooo, title, number, limitation)

# %% Reversion of the different lists
"Utility unknown"
aire.reverse()
gamma.reverse()
mach_function.reverse()
hotgas_temperature.reverse()
xcanauxre.reverse()
ycanauxre.reverse()
larg_canalre.reverse()
htre.reverse()
Areare.reverse()
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

# %% Preparation of the lists for CAO modelisation
"Changing the coordinates of the height of the canals (otherwise it is geometrically wrong)"
angles = []
newxhtre = []
newyhtre = []
for i in range(0, len(xcanauxre)):
    if i == 0:
        angle = 0
        angles.append(angle)
    elif i == (len(xcanauxre) - 1):
        angle = angles[i - 1]
        angles.append(angle)
    else:
        vect1 = (xcanauxre[i] - xcanauxre[i - 1]) / (
                (((ycanauxre[i] - ycanauxre[i - 1]) ** 2) + ((xcanauxre[i] - xcanauxre[i - 1]) ** 2)) ** 0.5)
        vect2 = (xcanauxre[i + 1] - xcanauxre[i]) / (
                (((ycanauxre[i + 1] - ycanauxre[i]) ** 2) + ((xcanauxre[i + 1] - xcanauxre[i]) ** 2)) ** 0.5)
        angle1 = degrees(acos(vect1))
        angle2 = degrees(acos(vect2))
        angle = angle2
        angles.append(angle)
    newx = xcanauxre[i] + htre[i] * sin(radians(angles[i]))
    newy = ycanauxre[i] + htre[i] * cos(radians(angles[i]))
    newxhtre.append(newx)
    newyhtre.append(newy)

"Checking the height of channels"
verification = []
for i in range(0, len(xcanauxre)):
    verifhtre = (((newxhtre[i] - xcanauxre[i]) ** 2) + ((newyhtre[i] - ycanauxre[i]) ** 2)) ** 0.5
    verification.append(verifhtre)
plt.plot(newxhtre, newyhtre)
plt.plot(xcanauxre, ycanauxre)
plt.title("Geometrical aspect of the canal")
plt.axis("equal")
plt.show()
plt.plot(xcanauxre, verification)
plt.title("Checking of the height of the generated canals")
plt.show()

# %% Writing the results of the study in a CSV file
"Writing the results in a CSV file"
Bar = ProgressBar(100, 30, "Writting results in CSV files    ")

file_name = "valuexport.csv"
file = open(file_name, "w")
writer = csv.writer(file)
writer.writerow(("Axe x moteur", "Diamètre moteur", "Aire gaz moteur", "Gamma gaz", "Nombre de mach", "Pression gaz",
                 "Pression totale", "Temperature gaz", "Axe x canaux", "Diamètre moteur+chemise", "Largeur canaux",
                 "Hauteur canaux", "Aire canaux", "Viscosité gaz", "Cp gaz", "Conductivité gaz", "Prandtl gaz",
                 "Coeff Hg", "Sigma", " Twg ", " Twl ", "Heat flux", "Tl", "Vitesse coolant", "Reynolds CH4",
                 "Coeff Hl", "Rho coolant", "Viscosité CH4", "Conductivité CH4", "Cp CH4", "Vitesse du coolant",
                 "Pression du coolant", "Conductivité de la paroi", "x hauteur réelle",
                 "y hauteur réelle"))
for i in range(0, len(xcanauxre)):
    writer.writerow(
        (x_value[i], y_value[i], aire_saved[i], gamma_saved[i], mach_function_saved[i], pressure_function[i],
         hotgas_temperature_saved[i], xcanauxre[i], ycanauxre[i], larg_canalre[i],
         htre[i], Areare[i], visc_function[i], cp_function[i], lamb_function[i], Prandtl_function[i],
         hg_function[i], Sig[i], inwall_temperature[i], outwall_temperature[i], fluxsolved[i], Tcoolant[i],
         Vitesse[i], Re_function[i], hlnormal[i], rho[i], visccoolant[i], condcoolant[i], Cpmeth[i],
         Vitesse[i], Pcoolant[i], LambdaTC[i], newxhtre[i], newyhtre[i]))
for i in range(len(xcanauxre), len(x_value)):
    writer.writerow(
        (x_value[i], y_value[i], aire_saved[i], gamma_saved[i], mach_function_saved[i], pressure_function[i],
         hotgas_temperature_saved[i], ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
         ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))
file.close()

# %% Writing the results of the study in a CSV file
"Writing the results of the study in a CSV file"
file_name = "geometry1.csv"
file = open(file_name, "w")
writer = csv.writer(file)
writer.writerow(("x hauteur réelle", "y hauteur réelle"))
for i in range(0, len(xcanauxre)):
    writer.writerow((newxhtre[i] * (-1000), newyhtre[i] * 1000))
for i in range(len(xcanauxre), len(x_value)):
    writer.writerow(
        (x_value[i], y_value[i], aire_saved[i], gamma_saved[i], mach_function_saved[i], pressure_function[i],
         hotgas_temperature_saved[i], ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
         ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))
file.close()

# %% Writing the results of the study in a CSV file
"Writing the results of the study in a CSV file"
file_name = "geometry2.csv"
file = open(file_name, "w")
writer = csv.writer(file)
writer.writerow(("Diamètre moteur+chemise", "x hauteur réelle"))
for i in range(0, len(xcanauxre)):
    writer.writerow((ycanauxre[i] * 1000, newxhtre[i] * (-1000)))
for i in range(len(xcanauxre), len(x_value)):
    writer.writerow(
        (x_value[i], y_value[i], aire_saved[i], gamma_saved[i], mach_function_saved[i], pressure_function[i],
         hotgas_temperature_saved[i], ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
         ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))
file.close()

Bar.update(100)
print()
print("█                                                                          █")
print("█                                                                          █")
print("███████████████████████████████████ FIN ████████████████████████████████████")

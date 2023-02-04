# %% CTE
"""
Created on Fri Nov 27 14:47:27 2020

@author: Julien

Rewritten: Mehdi

WARNING: This Python file was rewritten only for the Viserion_2023 project.
Any changes might affect the results.
"""
import csv
from pylab import *
from machsolve import Mach_solv
from pressuresolve import Pressure_solv
from temperaturesolve import Temperature_solv
from musolve import mu
from Tcorsolve import tempcorrige
from Canaux import canaux
from methan import *
from graphic3d import view3d
from heatequationsolve import *
from volume3d import *
from CoolProp.CoolProp import PropsSI
from IA import *
import matplotlib.pyplot as plt

print("██████████████████████████ Cool The Engine V 2.0.0 █████████████████████████")
print("█                  Innovative Propulsion Laboratory - IPL                  █")
# %% Engine initialisation
"Viserion settings"
plagex = "Viserion_X.txt"  # X coordinates of the Viserion
plagey = "Viserion_Y.txt"  # Y coordinates of the Viserion
plageinit = "Viserion_2023.txt"  # Viserion's parameters (found with CEA)

"Constant value"
lim22 = 600
size2 = 18
ange = -18.884
lim11 = 200
epso = 10
machtype = 0
limitation = 0.05

# %% Reading Viserion_2023.txt
"Reading of the file"
crinit = csv.reader(open(plageinit, "r"))
value = []
for row in crinit:
    value.append(row[1])

"Reading: Changeable variable according to CEA"
c_init = float(value[0])  # Sound velocity in the chamber
c_col = float(value[1])  # Sound velocity in the throat
debit_LOX = float(value[2])  # LOX debit
debit_LCH4 = float(value[3])  # CH4 debit
rho_init = float(value[4])  # Initial density of the gases
Pc = float(value[5])  # Pressure in the chamber
Tc = float(value[6])  # Combustion temperature (in the chamber?)
gamma_c = float(value[7])  # Gamma in the chamber
gamma_t = float(value[8])  # Gamma in the throat
gamma_e = float(value[9])  # Gamma at the exit
M = float(value[10])  # Molar mass of the gases
Cstar = float(value[11])  # caracteristic velocity

"Reading: Constant variable -- DO NOT CHANGE"
Dcol = float(value[12])  # Convergent radius of curvature
Rcol = float(value[13])  # Throat radius of curvature
Ac = float(value[14])  # Throat diameter
DiamCol = float(value[15])  # Throat diameter

# %% Import of the (X,Y) coordinates of the Viserion
"Reading the files"
crx = csv.reader(open(plagex, "r"))
cry = csv.reader(open(plagey, "r"))
x_value = []
y_value = []

print("█                                                                          █")
Bar = ProgressBar(100, 30, "Import of the coordinates       ")

"Importing X coordinates in a list"
for row in crx:
    a = float(row[0]) / 1000
    x_value.append(a)
ax = 100 / (len(x_value) - 1)
b = 0

"Importing Y coordinates in a list"
for row in cry:
    a = float(row[0]) / 1000
    y_value.append(a)
    b = b + ax
    Bar.update(b)

"Creation of the mesh"
R = []
for i in range(0, len(x_value) - 1, 1):
    maillage = abs(x_value[i] - x_value[i + 1])
    R.append(maillage)
R.append(R[-1])

# Plot of the upper profile of the engine
"""plt.figure(dpi=200)
plt.plot(x_value, y_value, color='black')
plt.title('Profile of the Viserion', color='black')
plt.show()"""

# Plot of the mesh density of the engine
"""
colooo = plt.cm.binary
inv = 1, 1, 1
view3d(inv, x_value, y_value, R, colooo, 'Mesh density', size2 - 2, limitation)
"""

print()

# %% Areas computation
"Computation of the cross-sectional areas of the engine"
long = len(x_value)
aire = []
Bar = ProgressBar(100, 30, "Computation of the areas        ")
aw = 100 / (long)
b = 0

for i in range(0, long, 1):
    a = pi * (y_value[i]) ** 2
    aire.append(a)
    b = b + aw
    Bar.update(b)

print()

# %% Adiabatic constant parametrization
"Computation of gamma  --  Not much information"
i = 0
a = 1
b = 1
while a == b:
    a = y_value[i]
    i = i + 1
    b = y_value[i]

gamma = []
for j in range(0, i, 1):
    gamma.append(gamma_c)

j = y_value.index(min(y_value))
k = j - i
c = gamma_c
a = -1
for m in range(0, k, 1):
    l = (gamma_c - gamma_t) / ((x_value[j] - x_value[i]) / abs(x_value[i + 1 + a] - x_value[i + a]))
    c = c - l
    a = a + 1
    gamma.append(c)

p = len(x_value) - j
c = gamma_t
a = -1
for q in range(0, p, 1):
    n = (gamma_t - gamma_e) / ((y_value[-1] - y_value[j]) / abs(y_value[j + 1 + a] - y_value[j + a]))
    c = c - n
    a = a + 1
    gamma.append(c)

# Plot of the gamma linearisation
"""
#print(gamma)
#print(len(gamma))
plt.figure(dpi=200)
plt.plot(x_value,gamma,color='gold')
plt.title("Gamma linearisation")
plt.show()
"""

# %% Mach number computation
"Computation of the initial velocity and mach number of the gases"
v_init = (debit_LOX + debit_LCH4) / (rho_init * aire[0])  # initial velocity of the gases
M_init = v_init / c_init
M1 = M_init
mach_function = [M_init]
b = 0
Bar = ProgressBar(100, 30, "Gamma, Mach number and pressure computations")
av = 100 / (long - 1)

"Mach number computations along the engine"
for i in range(0, long - 1, 1):
    A1 = aire[i]
    A2 = aire[i + 1]
    pos = i
    g = Mach_solv(A1, A2, M1, gamma[i], pos, machtype)
    mach_function.append(g)
    M1 = g
    b = b + av
    Bar.update(b)

# Plots of the Mach number in the engine (2D/3D)
plt.figure(dpi=200)
plt.plot(x_value, mach_function, color='gold')
plt.title("Mach number as a function of the engine axis")
plt.show()

colooo = plt.cm.Spectral
inv = 1, 1, 1
view3d(inv, x_value, y_value, mach_function, colooo, 'Mach number', size2 - 2, limitation)

print()

# %% Static pressure computation
"Static pressure computation"
pressure_function = []
pressure_function.append(Pc)
c = 0
Bar = ProgressBar(100, 30, "Static pressure computation")
ac = (long - 1) / 100

"Static pressure computations along the engine"
for i in range(0, long - 1, 1):
    if (i == long + 1):
        M1 = mach_function[i]
        M2 = mach_function[i]
    else:
        M1 = mach_function[i]
        M2 = mach_function[i + 1]
    P1 = pressure_function[i]
    P2 = Pressure_solv(M1, M2, P1, gamma[i])
    pressure_function.append(P2)
    c = c + ac
    Bar.update(c)

# Plot of the static pressure (2D/3D)
plt.figure(dpi=200)
plt.plot(x_value, pressure_function, color='gold')
plt.title("Static pressure as a function of the engine axis")
plt.show()

colooo = plt.cm.gist_rainbow_r
inv = 1, 1, 1
view3d(inv, x_value, y_value, pressure_function, colooo, 'Static pressure', size2 - 2, limitation)

print()

# %% Temperature computation
"Temperature computation"
temperature_function = []
temperature_function.append(Tc)
b = 0
Bar = ProgressBar(100, 30, "Temperature computation         ")
ay = 100 / (long - 1)

"Temperature computations along the engine"
for i in range(0, long - 1, 1):
    if (i == long + 1):
        M1 = mach_function[i]
        M2 = mach_function[i]
    else:
        M1 = mach_function[i]
        M2 = mach_function[i + 1]
    T1 = temperature_function[i]
    T2 = Temperature_solv(M1, M2, T1, gamma[i])
    temperature_function.append(T2)
    b = b + ay
    Bar.update(b)

"import temperature,gamma and mach number values in a list"
Tg_function = []
for i in range(0, long, 1):
    Tg = tempcorrige(temperature_function[i], gamma[i], mach_function[i])
    Tg_function.append(Tg)

# Plot of the temperature (2D/3D)
"""
plt.figure(dpi=200)
plt.plot(x_value,temperature_function,color='gold')
plt.title("Temperature as a function of the engine axis")
plt.show()

colooo = plt.cm.terrain_r
inv = 1, 1, 1
view3d(inv,x_value,y_value,temperature_function,colooo,'Temperature of the gases',size2-2,limitation)
"""

print()

# %% Total pressure computation
"Total pressure computation"
rho_function = []
rho_function.append(rho_init)
Celerite_function = [c_init]
Gaz_velocity = [v_init]
Ptotale_function = [pressure_function[0] + 0.5 * rho_init * v_init ** 2]
b = 0
Bar = ProgressBar(100, 30, "Total pressure computation      ")
ay = 100 / (long - 1)

"Total pressure computations along the engine"
for i in range(0, long - 1, 1):
    if (i == long + 1):
        M1 = mach_function[i]
        M2 = mach_function[i]
    else:
        M1 = mach_function[i]
        M2 = mach_function[i + 1]
    Cele = (gamma[i] * 434.47 * temperature_function[i]) ** 0.5
    Gaz_velo = Cele * mach_function[i]
    Celerite_function.append(Cele)
    Gaz_velocity.append(Gaz_velo)
    R1 = rho_function[i]
    R2 = Pressure_solv(M1, M2, R1, gamma[i])
    rho_function.append(R2)
    Ptot = pressure_function[i] + 0.5 * R2 * Gaz_velo ** 2
    Ptotale_function.append(Ptot)
    b = b + ay
    Bar.update(b)

# Plot of the total pressure (2D/3D)
"""
plt.plot(x_value,Ptotale_function,color='red')
plt.title("Pression totale des gaz en fonction de l'axe moteur")
plt.show()
"""

print()
print("█                                                                          █")

# %% Canal parameters
"""Number of canal and tore position"""
nbc = 40  # Number of canal
tore = 0.068  # Position of the tore from the throat (in m)

"Width of the canal"
lrg_c2 = 0.003  # Width of the canal in the high section of the chamber (in m)
lrg_c = 0.003  # Width of the canal in the low section of the chamber (in m)
lrg_col = 0.0035  # Width of the canal in the throat (in m)
lrg_div = 0.003  # Width of the canal in the divergent (in m)

"Height of the canal"
ht_c2 = 0.003  # Height of the canal in the high section of the chamber (in m)
ht_c = 0.003  # Height of the canal in the low section of the chamber (in m)
ht = 0.0035  # Height of the canal in the throat (in m)
ht_div = 0.003  # Height of the canal in the divergent (in m)

# %% Thickness
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
"Properties of the copper"
condcoeff1 = -0.065665283166
condcoeff2 = 421.82710859

# %% Coolant
"Properties of the CH4"
debit_LCH4 = debit_LCH4
fluid = "Methane"
rho_initCH4 = 425  # Volumic mass of the CH4 -- Not necessary
If_reg = debit_LCH4  # Total debit (in kg/s)
Tl_init = 111  # Initial temperature of the coolant (in K)
debit_total = If_reg / rho_initCH4  # Total volumic debit of the coolant (in m3/s)
Pl_init = 3700000  # Initial pressure of the coolant (in Pa)
Ro = 3  # UNKNOWN?????

# %% Computation
"""Methode 2"""
xcanauxre, ycanauxre, larg_canalre, Areare, htre, reste, epaiss_chemise = canaux(plagex, plagey, nbc, lrg_col, lrg_c,
                                                                                 lrg_div, ht, ht_c, ht_div, tore,
                                                                                 debit_total, n1, n2, n3, n4, e_col,
                                                                                 e_div, e_c, n5, n6, lrg_c2, ht_c2)
"""Methode 1"""
Bar = ProgressBar(100, 30, "Canal geometric computations    ")
Bar.update(100)
print()
print("█                                                                          █")
epaiss_chemise.reverse()
xcanauxre.reverse()
larg_canalre.reverse()
Areare.reverse()
htre.reverse()
ycanauxre.reverse()
fin = (len(xcanauxre))
Stemperature_function = temperature_function[:]
Saire = aire[:]
Smach_function = mach_function[:]
Sgamma = gamma[:]
while temperature_function.index(temperature_function[-1]) >= fin:
    temperature_function.pop()
    aire.pop()
    mach_function.pop()
    gamma.pop()

gamma.reverse()
mach_function.reverse()
aire.reverse()
temperature_function.reverse()


def mainsolver(Sig, b, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, ay, Pcoolant, LambdaTC, entropy):
    def downsolver(Sig, b, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, ay, Pcoolant, LambdaTC, entropy):
        pos = 0
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
        Burnout = []
        Burnout2 = []
        Celerite = []
        hlnormal = []
        error_D_ = []
        singpertes = [Pcoolant[0]]
        totalpha = 0
        Pcoolant2 = [Pcoolant[0]]
        Burnoutplus1 = []
        Burnoutplus2 = []
        phase = 0
        Quality = -1
        positioncol = ycanauxre.index(min(ycanauxre))
        for i in range(0, len(xcanauxre), 1):
            Lambda_tc = LambdaTC[i]
            x = xcanauxre[i]
            c = larg_canalre[i]
            Dhy = (2 * htre[i] * c) / (htre[i] + c)
            V = ((debit_LCH4 / (nbc * rho[i])) / Areare[i])
            Vitesse.append(V)
            Re = (V * Dhy * rho[i]) / visccoolant[i]
            Re_function.append(Re)
            Pr_cool = (visccoolant[i] * Cpmeth[i]) / condcoolant[i]
            T1 = temperature_function[i]
            viscosite, cp, lamb, Pr = mu(T1, M, gamma[i])
            visc_function.append(viscosite)
            cp_function.append(cp)
            lamb_function.append(lamb)
            Prandtl_function.append(Pr)
            A = aire[i]
            hg = ((0.026 / ((2 * ((Ac / pi) ** 0.5)) ** 0.2)) * (((viscosite ** 0.2) * cp) / (Pr ** 0.6)) * (
                    (Pc / Cstar) ** 0.8) * ((DiamCol / Dcol) ** 0.1) * ((Ac / A) ** 0.9)) * Sig[i]
            hg_function.append(hg)
            Tg = T1
            steff = 5.6697 * 10 ** (-8)
            emissivity = 0.02
            qr = emissivity * steff * (T1 ** 4)
            if (i + 1) == (len(xcanauxre)):
                EtalArea = (2 * c + 2 * htre[i]) * abs(xcanauxre[i] - xcanauxre[i - 1])
            else:
                EtalArea = (2 * c + 2 * htre[i]) * abs(xcanauxre[i + 1] - x)
            Gdeb = debit_LCH4 / (Areare[i] * nbc)
            if phase == 0 or phase == 2:
                f = (1.82 * log10(Re) - 1.64) ** (-2)
                Nu = ((f / 8) * (Re - 1000) * Pr_cool) / (1 + 12.7 * ((f / 8) ** 0.5) * ((Pr_cool ** (2 / 3)) - 1))
                hl = Nu * (condcoolant[i] / Dhy)
                Xqual = PropsSI("Q", "P", Pcoolant[i], "T", Tcoolant[i], fluid)
            else:
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
                HVP = PropsSI("H", "Q", 1, "P", Pcoolant[i], fluid)
                HLP = PropsSI("H", "Q", 0, "P", Pcoolant[i], fluid)
                HHV = HVP - HLP
                flux1 = flux + 0.5
                flux2 = flux
                while (abs(flux1 - flux2) / flux2) > 0.000001:
                    flux1 = (flux1 + flux2) / 2
                    Bo = 1000000 * flux1 / (Gdeb * HHV)
                    if Xqual < 0.6:
                        Nu = 12.46 * (Bo ** 0.544) * (We ** 0.035) * (Kp ** 0.614) * (Xeq ** 0.031)
                    else:
                        Nu = 0.00136 * (Bo ** (-1.442)) * (We ** 0.074)
                    hl = Nu * (condcoolant[i] / Dhy)
                    D = 2 * (ycanauxre[i] - epaiss_chemise[i])
                    d = (pi * (D + htre[i] + epaiss_chemise[i]) - nbc * c) / nbc
                    m = ((2 * hl) / (d * Lambda_tc)) ** 0.5
                    hl_cor = hl * ((nbc * c) / (pi * D)) + nbc * (
                            (2 * hl * Lambda_tc * (((pi * D) / nbc) - c)) ** 0.5) * ((tanh(m * htre[i])) / (pi * D))
                    hg = hg_function[i]
                    hl = hl_cor
                    Tl = Tcoolant[i]
                    e = epaiss_chemise[i]
                    L = Lambda_tc
                    mp.dps = 150
                    cx1 = Symbol('cx1')
                    cx2 = Symbol('cx2')
                    f1 = hg * (Tg - cx1) - (L / e) * (cx1 - cx2)
                    f2 = hl * (cx2 - Tl) - (L / e) * (cx1 - cx2)
                    x_, y_ = nsolve((f1, f2), (cx1, cx2), (900, 700))
                    flux2 = hl * (y_ - Tcoolant[i]) * 0.000001
            hl = Nu * (condcoolant[i] / Dhy)
            hlnormal.append(hl)
            D = 2 * (ycanauxre[i] - epaiss_chemise[i])
            d = (pi * (D + htre[i] + epaiss_chemise[i]) - nbc * c) / nbc
            m = ((2 * hl) / (d * Lambda_tc)) ** 0.5
            hl_cor = hl * ((nbc * c) / (pi * D)) + nbc * ((2 * hl * Lambda_tc * (((pi * D) / nbc) - c)) ** 0.5) * (
                    (tanh(m * htre[i])) / (pi * D))
            hlcor.append(hl_cor)
            hg = hg_function[i]
            hl = hlcor[i]
            Tl = Tcoolant[i]
            e = epaiss_chemise[i]
            L = Lambda_tc
            mp.dps = 150
            cx1 = Symbol('cx1')
            cx2 = Symbol('cx2')
            f1 = hg * (Tg - cx1) - (L / e) * (cx1 - cx2)
            f2 = hl * (cx2 - Tl) - (L / e) * (cx1 - cx2)
            x_, y_ = nsolve((f1, f2), (cx1, cx2), (900, 700))
            inwall_temperature.append(x_)
            outwall_temperature.append(y_)
            flux = hl * (y_ - Tcoolant[i]) * 0.000001
            fluxsolved.append(flux)
            Tw = inwall_temperature[i]
            Ts = temperature_function[positioncol]
            Mak = mach_function[i]
            sigm = 1 / ((((Tw / (2 * Ts)) * (1 + (((gamma[i] - 1) / 2) * (Mak ** 2))) + 0.5) ** (0.68)) * (
                    (1 + (((gamma[i] - 1) / 2) * (Mak ** 2))) ** 0.12))
            Sig.append(sigm)
            Lambdatc = (condcoeff1 * ((Tw + outwall_temperature[i]) * 0.5) + condcoeff2)  # 16.87
            LambdaTC.append(Lambdatc)
            if (i == xcanauxre.index(xcanauxre[-1])):
                Distance = ((xcanauxre[i - 1] - xcanauxre[i]) ** 2 + (ycanauxre[i - 1] - ycanauxre[i]) ** 2) ** 0.5
                xa = Distance
                ya = (2 * pi * ycanauxre[i - 1]) / nbc
                za = (2 * pi * ycanauxre[i]) / nbc
                perim = (2 * pi * ((A / pi) ** 0.5)) / nbc
                dA = xa * (2 * c + 2 * htre[i])
                timer = (((xcanauxre[i - 1] - xcanauxre[i]) ** 2 + (ycanauxre[i - 1] - ycanauxre[i]) ** 2) ** 0.5) / V
            else:
                Distance = ((xcanauxre[i + 1] - xcanauxre[i]) ** 2 + (ycanauxre[i + 1] - ycanauxre[i]) ** 2) ** 0.5
                xa = Distance
                ya = (2 * pi * ycanauxre[i + 1]) / nbc
                za = (2 * pi * ycanauxre[i]) / nbc
                perim = (2 * pi * ((A / pi) ** 0.5)) / nbc
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
            if phase == 0:
                H1 = PropsSI("H", "P", Pcoolant[i], "T", Tcoolant[i], fluid)
                H2 = H1 + Q / debitmass
                Quality = PropsSI("Q", "H", H2, "P", Pcoolant[i + 1], fluid)
                if Quality > 0 and Quality < 1:
                    phase = 1
                else:
                    phase = 0
            else:
                H1 = H2
                H2 = H1 + Q / debitmass
                Quality = PropsSI("Q", "H", H2, "P", Pcoolant[i + 1], fluid)
                if Quality > 0 and Quality < 1:
                    phase = 1
                else:
                    phase = 2
            if phase == 0 or phase == 2:
                Tcoolant.append(Tfu)
            else:
                Tcoolant.append(Tcoolant[i])
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
            PMPa = Pcoolant[i] / 1000000
            if PMPa <= 5.5:
                first = 0.82281553398 * PMPa + 0.91213592233
            else:
                first = 0.076974683097 * PMPa + 4.2238162494
            newBurnout = first + 0.334433962 * (V - 10)
            Burnout.append(newBurnout)
            if phase == 0 or phase == 1:
                DTemp = float(DeltaT(Pcoolant[i], Tcoolant[i], fluid))
                Burn2x = (V * 3.281) * (DTemp * (9 / 5))
                Burn2 = (4.68521 * 10 ** (-3)) * Burn2x + 1.00130 * 10
                Burn3 = 16.270 * (0.2598 + 0.0004134 * DTemp * (9 / 5) * ((V * 3.281) ** 0.9))
            else:
                Burn2 = 0
                Burn3 = 0
            Burnoutplus1.append(Burn2)
            Burnoutplus2.append(Burn3)
            b = b + ay
            Bar.update(b)
        Burnout2.append(Burnoutplus1)
        Burnout2.append(Burnoutplus2)
        return hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, inwall_temperature, \
               outwall_temperature, fluxsolved, Sig, b, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, \
               Vitesse, Pcoolant, LambdaTC, Burnout, Burnout2, Celerite, hlnormal, error_D_, singpertes, Pcoolant2

    hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, inwall_temperature, \
    outwall_temperature, fluxsolved, Sig, b, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, Vitesse, \
    Pcoolant, LambdaTC, Burnout, Burnout2, Celerite, hlnormal, error_D_, singpertes, Pcoolant2 = downsolver(
        Sig, b, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, ay, Pcoolant, LambdaTC, entropy)
    return hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, inwall_temperature, \
           outwall_temperature, fluxsolved, Sig, b, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, \
           Vitesse, Pcoolant, LambdaTC, Burnout, Burnout2, Celerite, hlnormal, error_D_, singpertes, Pcoolant2


Sig = []
Tcoolant = []
Pcoolant = []
visccoolant = []
condcoolant = []
Cpmeth = []
rho = []
LambdaTC = []
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
Bar = ProgressBar(100, 30, "Résolution globale en cours...  ")
ay = 100 / ((1 + rep) * len(xcanauxre))
hlcor, visc_function, cp_function, lamb_function, Prandtl_function, hg_function, inwall_temperature, \
outwall_temperature, fluxsolved, Sig, b, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, Vitesse, \
Pcoolant, LambdaTC, Burnout, Burnout2, Celerite, hlnormal, error_D_, singpertes, Pcoolant2 = mainsolver(
    Sig, b, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, ay, Pcoolant, LambdaTC, entropy)
for i in range(0, rep, 1):
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
    outwall_temperature, fluxsolved, Sig, b, Re_function, Tcoolant, visccoolant, condcoolant, Cpmeth, rho, Vitesse, \
    Pcoolant, LambdaTC, Burnout, Burnout2, Celerite, hlnormal, error_D_, singpertes, Pcoolant2 = mainsolver(
        Sig, b, rho, Tcoolant, visccoolant, condcoolant, Cpmeth, ay, Pcoolant, LambdaTC, entropy)

print()
print("█                                                                          █")

# %% Display of the first results
"Display of the results"
colooo = plt.cm.magma
inv = 0, 0, 0
view3d(inv, xcanauxre, ycanauxre, inwall_temperature, colooo, "Wall temperature on the gas side", size2, limitation)
Cel03 = []
for x in Celerite:
    x = x * 0.3
    Cel03.append(x)

"""plt.figure(dpi=200)
plt.plot(xcanauxre, Re_function, color='blue')
plt.title("Reynolds number as a function of the engine axis")
plt.show()
plt.figure(dpi=200)
plt.plot(xcanauxre, hlcor, color='blue', label='Hl corrigé')
plt.plot(xcanauxre, hlnormal, color='cyan', label='Hl')
plt.title("Convection coefficient Hl as a function of the engine axis")
plt.legend()
plt.show()
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
plt.show()"""

plt.figure(dpi=200)
plt.plot(xcanauxre, inwall_temperature, color='magenta', label='Twg')
plt.plot(xcanauxre, outwall_temperature, color='orangered', label='Twl')
plt.title('Wall temperature')
plt.legend()
plt.show()
plt.figure(dpi=200)
plt.plot(xcanauxre, fluxsolved, color='darkviolet', label='Heat flux')
plt.plot(xcanauxre, Burnout2[0], color='black', label='Burn out NASA')
plt.plot(xcanauxre, Burnout2[1], color='red', label='Burn out Pratt & Witney')
plt.title('Burn out')
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
view3d(inv, xcanauxre, ycanauxre, fluxsolved, colooo, "Burnout 3D", size2, limitation)
colooo = plt.cm.coolwarm
inv = 0, 0, 0
view3d(inv, xcanauxre, ycanauxre, Tcoolant, colooo, "Temperature of the coolant", size2, limitation)

# %% Flux computation in 2D and 3D
"Computation for 2D graph"
reste.reverse()
pas = reste[-1] + larg_canalre[-1]
epaisseur = e_c
hauteur = htre[-1]
largeur = larg_canalre[-1]
Hg = hg_function[-1]
Tg = temperature_function[-1]
Hl = hlnormal[-1]
Tl = Tcoolant[-1]
dx = 0.00004  # *3.5
lamb = LambdaTC[-1]
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 5, 1, 1)
poscol = ycanauxre.index(min(ycanauxre))
pas = reste[poscol] + larg_canalre[poscol]
epaisseur = e_col
hauteur = htre[poscol]
largeur = larg_canalre[poscol]
Hg = hg_function[poscol]
Tg = temperature_function[poscol]
Hl = hlnormal[poscol]
Tl = Tcoolant[poscol]
dx = 0.000025  # *3.5
lamb = LambdaTC[poscol]
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 15, 1, 2)
pas = reste[0] + larg_canalre[0]
epaisseur = e_div
hauteur = htre[0]
largeur = larg_canalre[0]
Hg = hg_function[0]
Tg = temperature_function[0]
Hl = hlnormal[0]
Tl = Tcoolant[0]
dx = 0.00004
lamb = LambdaTC[0]
t3d = carto2D(pas, epaisseur, hauteur, largeur, dx, Hg, lamb, Tg, Hl, Tl, 5, 1, 1)

"Computation for 3D graph"
eachT = []
for i in range(0, len(xcanauxre), 1):
    # print(i)
    lim1 = 0
    lim2 = 650
    pas = reste[i] + larg_canalre[i]
    epaisseur = epaiss_chemise[i]
    hauteur = htre[i]
    largeur = larg_canalre[i]
    Hg = hg_function[i]
    Tg = temperature_function[i]
    Hl = hlnormal[i]
    Tl = Tcoolant[i]
    if i <= lim2 and i >= lim1:
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
temperature_function.reverse()
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

# %% Preparation of the lists for CAO
"Changing the coordinates of the height of the canals (otherwise it is geometrically wrong)"
angles = []
newxhtre = []
newyhtre = []
for i in range(0, len(xcanauxre), 1):
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

"checking the height"
verification = []
for i in range(0, len(xcanauxre), 1):
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
for i in range(0, len(xcanauxre), 1):
    writer.writerow((x_value[i], y_value[i], Saire[i], Sgamma[i], Smach_function[i], pressure_function[i],
                     Ptotale_function[i], Stemperature_function[i], xcanauxre[i], ycanauxre[i], larg_canalre[i],
                     htre[i], Areare[i], visc_function[i], cp_function[i], lamb_function[i], Prandtl_function[i],
                     hg_function[i], Sig[i], inwall_temperature[i], outwall_temperature[i], fluxsolved[i], Tcoolant[i],
                     Vitesse[i], Re_function[i], hlnormal[i], rho[i], visccoolant[i], condcoolant[i], Cpmeth[i],
                     Vitesse[i], Pcoolant[i], LambdaTC[i], newxhtre[i], newyhtre[i]))
for i in range(len(xcanauxre), len(x_value), 1):
    writer.writerow((x_value[i], y_value[i], Saire[i], Sgamma[i], Smach_function[i], pressure_function[i],
                     Ptotale_function[i], Stemperature_function[i], ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                     ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))
file.close()

# %% Writing the results of the study in a CSV file
"Writing the results of the study in a CSV file"
file_name = "geometry1.csv"
file = open(file_name, "w")
writer = csv.writer(file)
writer.writerow(("x hauteur réelle", "y hauteur réelle"))
for i in range(0, len(xcanauxre), 1):
    writer.writerow((newxhtre[i] * (-1000), newyhtre[i] * 1000))
for i in range(len(xcanauxre), len(x_value), 1):
    writer.writerow((x_value[i], y_value[i], Saire[i], Sgamma[i], Smach_function[i], pressure_function[i],
                     Ptotale_function[i], Stemperature_function[i], ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                     ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))
file.close()

# %% Writing the results of the study in a CSV file
"Writing the results of the study in a CSV file"
file_name = "geometry2.csv"
file = open(file_name, "w")
writer = csv.writer(file)
writer.writerow(("Diamètre moteur+chemise", "x hauteur réelle"))
for i in range(0, len(xcanauxre), 1):
    writer.writerow((ycanauxre[i] * 1000, newxhtre[i] * (-1000)))
for i in range(len(xcanauxre), len(x_value), 1):
    writer.writerow((x_value[i], y_value[i], Saire[i], Sgamma[i], Smach_function[i], pressure_function[i],
                     Ptotale_function[i], Stemperature_function[i], ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                     ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))
file.close()

print()
print("█                                                                          █")
print("███████████████████████████████████ FIN ████████████████████████████████████")

from CoolProp.CoolProp import PropsSI
import numpy as np


def rhoCH4(P, T, fluid):
    rho = PropsSI("D", "T", T, "P", P, fluid)
    return rho


def cpCH4(P, T, fluid):
    Cp = PropsSI("C", "P", P, "T", T, fluid)
    return Cp


def condCH4(P, T, fluid):
    conductivity = PropsSI("L", "T", T, "P", P, fluid)
    return conductivity


def viscCH4(P, T, fluid):
    viscosity = PropsSI("V", "T", T, "P", P, fluid)
    return viscosity


def sonCH4(P, T, fluid):
    celerite = PropsSI("A", "T", T, "P", P, fluid)
    return celerite


def entCH4(P, T, fluid):
    ent = PropsSI("S", "T", T, "P", P, fluid)
    return ent


def pressureCH4(S, T, fluid):
    h1 = PropsSI('P', 'T', T, 'S', S, fluid)
    return h1


def DeltaT(P, T, fluid):
    temp = PropsSI("T", "P", P, "Q", 0, fluid)
    Delta_T = temp - T
    return Delta_T


def pvk(Re, lambi):
    laaam = []
    errooor = []
    for i in range(0, 10000, 1):
        inn = 1 / 10000000
        lambi = lambi - inn
        error = (1 / (lambi ** 0.5)) - 2 * np.log10(Re * (lambi ** 0.5)) + 0.8
        laaam.append(lambi)
        errooor.append(abs(error))
    for i in range(0, 10000, 1):
        inn = 1 / 10000000
        lambi = lambi + inn
        error = (1 / (lambi ** 0.5)) - 2 * np.log10(Re * (lambi ** 0.5)) + 0.8
        laaam.append(lambi)
        errooor.append(abs(error))

    solute = min(errooor)
    isolute = errooor.index(solute)
    solution = laaam[isolute]
    return solution, solute

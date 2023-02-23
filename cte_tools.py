from sympy import Symbol, nsolve
import sympy as mp
import numpy as np
import scipy.optimize as opt

def mach_solv(area_1, area_2, mach_1, gamma, position, machtype):
    def solveur(area_1, area_2, mach_1, gamma, ome):
        """
        if 320 < position < 370:
            a = 10
        else:
            a = 1
        """
        if area_1 == area_2:
            solution_mach = mach_1
        else:
            part_2 = (area_1 / area_2) * (mach_1 / ((1 + ((gamma - 1) / 2) * mach_1 * mach_1) ** ome))
            M2 = mach_1
            liste = []
            mach = []
            for i in range(0, 2000):
                M2 += 0.00001
                part_1 = M2 * ((1 + (((gamma - 1) / 2) * M2 * M2)) ** (-ome))
                liste.append(abs(part_1 - part_2))
                mach.append(M2)
                # print("For mach =",M2,"diff=",diff)

            solution_mach = mach[liste.index(min(liste))]
        # print("Le nombre de mach est ",solution_mach,"avec une erreur de l'ordre de", min(liste))
        return solution_mach

    #  er=(1/0.07)*part_2**(ome/(2*ome+1))
    #  er2=part_2**ome
    #  print(er)
    #  print(er2)
    def solving(area_1, area_2, mach_1, gamma, ome):
        mp.dps = 1
        cx = Symbol('cx')
        f1 = (area_1 / area_2) * (mach_1 / ((1 + ((gamma - 1) / 2) * mach_1 * mach_1) ** ome)) - cx * (
                (1 + (((gamma - 1) / 2) * cx * cx)) ** (-ome))
        M2 = nsolve(f1, cx, mach_1, verify=False)
        return M2

    def solving_bis(area_1, area_2, mach_1, gamma, ome):
        def f1(x):
            f = (area_1 / area_2) * (mach_1 / ((1 + ((gamma - 1) / 2) * mach_1 * mach_1) ** ome)) - x * (
                (1 + (((gamma - 1) / 2) * x * x)) ** (-ome))
            return f
        M2 = opt.fsolve(func=f1, x0=mach_1)
        return M2

    ome = (gamma + 1) / (2 * (gamma - 1))
    if machtype == 0:
        final = solveur(area_1, area_2, mach_1, gamma, ome)
    elif(machtype == 1):
        final = solving(area_1, area_2, mach_1, gamma, ome)
    else:
        final = solving_bis(area_1, area_2, mach_1, gamma, ome)
    return final


def pressure_solv(M1, M2, P1, gamma):
    """
    Compute hot gas pressure at next point, given mach numbers and previous pressure
    """

    part1 = (gamma / (gamma - 1)) * np.log(M1 * M1 + (2 / (gamma - 1)))
    part2 = (gamma / (gamma - 1)) * np.log(M2 * M2 + (2 / (gamma - 1)))
    part3 = np.log(P1)

    return np.exp(part1 - part2 + part3)


def temperature_solv(M1, M2, T1, gamma):
    """
    Compute hot gas temperature at next point, given mach numbers and previous temperature
    """

    part1 = -np.log(abs(((gamma - 1) * M1 * M1) + 2))
    part2 = -np.log(abs(((gamma - 1) * M2 * M2) + 2))
    part3 = np.log(T1)

    return np.exp(-part1 + part2 + part3)


def tempcorrige(temp_original, gamma, mach_number):
    """
    Correct the hot gas temperature using a correlation from [INSERT SOURCE]
    """

    Pr = 4 * gamma / (9 * gamma - 5)
    temp_corrected = temp_original * (1 + (Pr ** 0.33) * (
            (gamma - 1) / gamma) * mach_number ** 2) / (1 + (
            (gamma - 1) / gamma) * mach_number ** 2)

    return temp_corrected


def conductivity(Twg: float, Twl: float, material_name: str):
    """
    Compute the conductivity of the chamber wall, given temperature and material
    """

    T_avg = (Twg + Twl) * 0.5
    if material_name == "pure copper":
        return -0.065665 * T_avg + 421.82
    if material_name == "cucrzr":
        return -0.0269 * T_avg + 365.33
    if material_name == "inconel":
        return 0.0138 * T_avg + 5.577


def hotgas_properties(gas_temp, molar_mass_, gamma):
    """
    Computes the properties of the hot gases according to [INSERT SOURCE].
    """

    dyn_viscosity = 17.858 * (46.61 * 10 ** (-10)) * (molar_mass_ ** 0.5) * ((9 / 5) * gas_temp) ** 0.6
    Cp = 8314 * gamma / ((gamma - 1) * molar_mass_)
    Lamb = dyn_viscosity * (Cp + (8314.5 / (4 * molar_mass_)))
    Pr = 4 * gamma / (9 * gamma - 5)

    return dyn_viscosity, Cp, Lamb, Pr


def flux_equations(guess, *data):
    """
    Used by scipy.optimize.fsolve() to compute hot and cold wall temperature.
    """

    t_hot, t_cold = guess  # Initial guess
    hg, hl, t_g, t_c, wall_conductivity, wall_thickness = data

    # System of equations to solve
    f1 = hg * (t_g - t_hot) - (wall_conductivity / wall_thickness) * (t_hot - t_cold)
    f2 = hl * (t_cold - t_c) - (wall_conductivity / wall_thickness) * (t_hot - t_cold)

    return [f1, f2]


def darcy_weisbach(Dhy, Re, roughness):
    friction_factor_1 = 1e-3
    friction_factor_2 = (1 / (-2 * np.log10(
        ((roughness / (Dhy * 3.7)) + 2.51 / (Re * (friction_factor_1 ** 0.5)))))) ** 2

    while abs((friction_factor_1 / friction_factor_2) - 1) > 0.0000001:
        friction_factor_1 = friction_factor_2
        friction_factor_2 = (1 / (-2 * np.log10(
            ((roughness / (Dhy * 3.7)) + 2.51 / (Re * (friction_factor_1 ** 0.5)))))) ** 2

    return friction_factor_2

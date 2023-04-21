import numpy as np


def mach_solv(area_1, area_2, mach_1, gamma):
    if area_1 == area_2:
        solution_mach = mach_1
    else:
        ome = (gamma + 1) / (2 * (gamma - 1))
        part_2 = (area_1 / area_2) * (mach_1 / ((1 + ((gamma - 1) / 2) * mach_1 * mach_1) ** ome))
        mach_2 = mach_1
        liste = []
        mach = []
        # Search of the mach_2 for which part_1 is minimum (750 iterations was ideal when tested)
        for i in range(0, 750):
            mach_2 += 0.00001
            part_1 = mach_2 * ((1 + (((gamma - 1) / 2) * mach_2 * mach_2)) ** (-ome))
            liste.append(abs(part_1 - part_2))
            mach.append(mach_2)
        solution_mach = mach[liste.index(min(liste))]
    return solution_mach


def pressure_solv(mach_1, mach_2, pressure_1, gamma):
    """
    Compute hot gas pressure at next point, given mach numbers and previous pressure
    """

    part1 = (gamma / (gamma - 1)) * np.log(mach_1 * mach_1 + (2 / (gamma - 1)))
    part2 = (gamma / (gamma - 1)) * np.log(mach_2 * mach_2 + (2 / (gamma - 1)))
    part3 = np.log(pressure_1)

    return np.exp(part1 - part2 + part3)


def temperature_hotgas_solv(mach_1, mach_2, temperature_1, gamma):
    """
    Compute hot gas temperature at next point, given mach numbers and previous temperature
    """

    part1 = np.log(abs(((gamma - 1) * mach_1 * mach_1) + 2))
    part2 = np.log(abs(((gamma - 1) * mach_2 * mach_2) + 2))
    part3 = np.log(temperature_1)

    return np.exp(part1 - part2 + part3)


def tempcorrige_pempie(temp_original, gamma, mach):
    """
    Compute the recovery temperature (adiabatic wall temperature) [P. Pempie]
    """

    Pr = 4 * gamma / (9 * gamma - 5)
    recovery_temp = temp_original * (
            (1 + (Pr ** 0.33) * ((gamma - 1) / 2) * mach ** 2) / (1 + ((gamma - 1) / 2) * mach ** 2))

    return recovery_temp


def conductivity(Twg: float, Twl: float, material_name: str):
    """
    Compute the conductivity of the chamber wall, given temperature and material
    """

    T_avg = (Twg + Twl) / 2
    if material_name == "pure copper":
        return -0.065665 * T_avg + 421.82
    if material_name == "cucrzr":
        return -0.0269 * T_avg + 365.33
    if material_name == "inconel":
        return 0.0138 * T_avg + 5.577


def hotgas_properties(gas_temp, molar_mass, gamma):
    """
    Computes the properties of the hot gases according to [INSERT SOURCE].
    """

    dyn_viscosity = 17.858 * (46.61 * 10 ** (-10)) * (molar_mass ** 0.5) * ((9 / 5) * gas_temp) ** 0.6
    Cp = 8314 * gamma / ((gamma - 1) * molar_mass)
    Lamb = dyn_viscosity * (Cp + (8314.5 / (4 * molar_mass)))
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

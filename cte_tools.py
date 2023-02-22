from sympy import Symbol, nsolve
import sympy as mp
import numpy as np


def mach_solv(A1, A2, M1, gamma, pos, machtype):
    def solveur(A1, A2, M1, gamma, ome):
        def part2(A1, A2, M1, gamma, ome):
            part2 = (A1 / A2) * (M1 / ((1 + ((gamma - 1) / 2) * M1 * M1) ** ome))

            return part2

        def part1(M2, ome, gamma):
            part1 = M2 * ((1 + (((gamma - 1) / 2) * M2 * M2)) ** (-ome))

            return part1

        """
        def mini(liste):
            mini = liste[0][0]
            for i in liste:
                if i[0] <= mini:
                    mini = i[0]
            return mini
        """
        if 320 < pos < 370:
            a = 10
        else:
            a = 1
        if A1 == A2:
            solution_mach = M1
        else:
            part_2 = part2(A1, A2, M1, gamma, ome)
            M2 = M1
            liste = []
            mach = []
            for i in range(0, 2000, 1):
                M2 += 0.00001
                part_1 = part1(M2, ome, gamma)
                diff = abs(part_1 - part_2)
                liste.append(diff)
                mach.append(M2)
                # print("For mach =",M2,"diff=",diff)

            soluce = min(liste)
            indice = liste.index(soluce)
            solution_mach = mach[indice]
        # print("Le nombre de mach est ",solution_mach,"avec une erreur de l'ordre de",soluce)
        return solution_mach

    #  er=(1/0.07)*part_2**(ome/(2*ome+1))
    #  er2=part_2**ome
    #  print(er)
    #  print(er2)
    def solving(A1, A2, M1, gamma, ome):
        mp.dps = 1
        cx = Symbol('cx')
        f1 = (A1 / A2) * (M1 / ((1 + ((gamma - 1) / 2) * M1 * M1) ** ome)) - cx * (
                (1 + (((gamma - 1) / 2) * cx * cx)) ** (-ome))
        M2 = nsolve(f1, cx, M1, verify=False)
        # print(M2)
        return M2

    """
    A1=1
    A2=3.7402
    M1=1
    gamma=1.1494
    """
    ome = (gamma + 1) / (2 * (gamma - 1))
    # final = solving(A1,A2,M1,gamma,ome)
    if machtype == 0:
        final = solveur(A1, A2, M1, gamma, ome)
    else:
        final = solving(A1, A2, M1, gamma, ome)
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
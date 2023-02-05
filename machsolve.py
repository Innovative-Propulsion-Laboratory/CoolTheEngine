# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 12:23:55 2020

@author: julie
"""
from sympy import Symbol, nsolve, Eq
import sympy as mp
from math import *


def Mach_solv(A1, A2, M1, gamma, pos, machtype):
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

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 15:08:38 2020

@author: julie
"""
import numpy as np


def primitive(x, gamma):
    a = (gamma / (gamma - 1)) * np.log(x * x + (2 / (gamma - 1)))
    return a


def pressure_solv(M1, M2, P1, gamma):
    part1 = primitive(M1, gamma)
    part2 = primitive(M2, gamma)
    part3 = np.log(P1)
    somme = part1 - part2 + part3
    P2 = np.exp(somme)
    return P2

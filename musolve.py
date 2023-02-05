# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:02:09 2020

@author: julie
"""


def mu(Tg, M, gamma):
    mu = 17.858 * (46.61 * 10 ** (-10)) * (M ** 0.5) * ((9 / 5) * Tg) ** 0.6
    Cp = 8314 * gamma / ((gamma - 1) * M)
    Lamb = mu * (Cp + (8314.5 / (4 * M)))
    Pr = 4 * gamma / (9 * gamma - 5)

    return mu, Cp, Lamb, Pr

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 21:29:44 2020

@author: julie
"""


def tempcorrige(temp_original, gamma, mach_number):
    Pr = 4 * gamma / (9 * gamma - 5)
    temp_corrected = temp_original * (1 + (Pr ** 0.33) * (
            (gamma - 1) / gamma) * mach_number ** 2) / (1 + (
            (gamma - 1) / gamma) * mach_number ** 2)
    return temp_corrected

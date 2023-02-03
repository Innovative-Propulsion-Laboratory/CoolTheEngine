# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 16:41:32 2020

@author: julien
"""

from math import *

def primitive(x,gamma):
    a=-log(abs(((gamma-1)*x*x)+2))
    return a
def Temperature_solv(M1,M2,T1,gamma):

    part1=primitive(M1,gamma)
    part2=primitive(M2,gamma)
    part3=log(T1)
    somme=-part1+part2+part3
    T2=exp(somme)
    return T2

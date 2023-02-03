# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 15:08:38 2020

@author: julie
"""
from math import *

def primitive(x,gamma):
    a=(gamma/(gamma-1))*log(x*x+(2/(gamma-1)))
    return a
def Pressure_solv(M1,M2,P1,gamma):
    part1=primitive(M1,gamma)
    part2=primitive(M2,gamma)
    part3=log(P1)
    somme=part1-part2+part3
    P2=exp(somme)
    return P2
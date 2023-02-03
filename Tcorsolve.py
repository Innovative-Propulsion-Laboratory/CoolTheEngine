# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 21:29:44 2020

@author: julie
"""

def tempcorrige(Tc,gamma,M):
    Pr=4*gamma/(9*gamma-5)
    Tg=Tc*(1+(Pr**0.33)*((gamma-1)/gamma)*M**2)/((1+((gamma-1)/gamma)*M**2))
    return Tg
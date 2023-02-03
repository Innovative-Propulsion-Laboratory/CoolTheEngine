# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 20:58:43 2021

@author: Luke 
"""

from __future__ import print_function
from CoolProp import AbstractState
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
import CoolProp.CoolProp as CoolProp
from CoolProp.HumidAirProp import HAPropsSI
"""
def viscCH4(P, T):
    #a = -6.060138165459480*10**(-11)*(T**6) + 5.085749710499200*10**(-8)*(T**5)  - 1.758382765515110*10**(-5)*(T**4)  + 3.204601949469970*10**(-3)*(T**3)  - 3.245501191668040*10**(-1)*(T**2)  + 1.731043198038020*10**(1)*(T**1)  - 3.795137193689100*10**2
    a= 1.529997958703050*10**(4)*T**(-2.423841635545800)
    b = -3.140970042392020*10**(-4)*(T**3) + 1.445484728370650*10**(-1)*(T**2) - 2.321753140731190*10**(1)*(T**1) + 1.345694014083140*10**3
    viscosity = (a*P + b)/1000000
    return viscosity

def rhoCH4(P, T):
    a = 6.335195878714910000000000000000*10**(-13)*(T**4) - 3.535795199395190000000000000000*10**(-10)*(T**3)+ 7.257389979897370000000000000000*10**(-8)*(T**2) - 6.486604886852970000000000000000*10**(-6)*T  + 2.127015847172960000000000000000*10**(-4)
    b = y = -2.337950803093090000000000000000*10**(-11)*(T**4) + 1.598526276369130000000000000000*10**(-8)*(T**3) - 3.793564868296250000000000000000*10**(-6)*(T**2)  + 3.774283292900870000000000000000*10**(-4)*T - 1.337860824899320000000000000000*10**(-2)
    c = y = -1.563339539438240000000000000000*10**(-8)*(T**4) + 8.034106039224920000000000000000*10**(-6)*(T**3) - 1.529616606992710000000000000000*10**(-3)*(T**2)  + 1.279049749337270000000000000000*10**(-1)*T - 3.968076288567900000000000000000

    d = y = 1.867913752472360000000000000000*10**(-6)*(T**4) - 9.844394864755200000000000000000*10**(-4)*(T**3) + 1.920987068785240000000000000000*10**(-1)*(T**2)  - 1.644112792794540000000000000000*10**(1)*T + 5.208016488611150000000000000000*10**(2)

    e = y = -5.805421368750350000000000000000*10**(-5)*(T**4) + 3.059318727480490000000000000000*10**(-2)*(T**3)- 5.967158147197240000000000000000*10**(-0)*(T**2)  + 5.085390694165220000000000000000*10**(2)*T- 1.551964247941030000000000000000*10**(4)

    rho = a*P**4 + b*P**3 + c*P**2 + d*P + e
    return rho

def condCH4(P, T):
    P=P/10
    a = (-7.538817638395410*10**(-7))*(T**3) + 3.054405390767840*10**(-4)*(T**2) - 4.068258623329730*10**(-2)*T + 1.782489859018910
    b = (1.5873984590*10**(-5))*(T**3) - 6.3850885518*10**(-3)*(T**2) + 8.4474478562*10**(-1)*T- 3.6780455766*10
    c =(-1.2135423922*10**(-4))*(T**3) + 4.8397105054*10**(-2)*(T**2) - 6.3535295023*T + 2.7470589687*100
    d =  (4.1316357026*10**(-4))*(T**3) - 1.6333362451*10**(-1)*(T**2) + 2.1296843774*10**(1)*(T) - 9.1487888381*100

    e = (-5.1364955031*10**(-4))*(T**3) + 1.9896992462*10**(-1)*(T**2) - 2.6844576132*10**(1)*T + 1.4130299162*1000
    Conductivity = (a*(P**4) + b*(P**3) + c*(P**2) + d*P + e)/1000
    return Conductivity

def cpCH4(T):
    #Cp =  1.116328788694580*10**(-12)*T**5 - 2.293819634094910*10**(-9)*T**4 +1.614578578737360*10**(-6)*T**3 - 3.878560932626700*10**(-4)*T**2 + 3.606597424639760*10**(-2)*T + 3.216651147426510*10
    #Cp=Cp/0.01604
    Cp = 7.2632015700*10**(-4)*T**4 - 3.9387427351*10**(-1)*T**3 + 7.9989764726*10**(1)*T**2 - 7.1916531769*10**(3)*T + 2.4464428487E+05
    return Cp
#y = E-12x5 -  1,614578578737360E-06x3 - 3,878560932626700E-04x2 + 3,606597424639760E-02x + 3,216651147426510E+01
"""
"""
P = 35 #atm
T = 111 #K
Conductivity = condCH4(P, T)
Viscosity = viscCH4(P, T)
Density = rhoCH4(P, T)
Cp=cpCH4(T)
print(Conductivity, Viscosity, Density,Cp)


"---------------------------------------------------------------------------"

#y= -7,539778739778620E-07x3 + 3,055328523328470E-04x2 - 4,070232130832050E-02x + 1,783652380952340E+00
#y = 1,5871899952E-05x3 - 6,3841891294E-03x2 + 8,4461739490E-01x - 3,6774523810E+01
#y = -5,1390630111E-04x3 + 1,9908271669E-01x2 - 2,6860841251E+01x + 1,4138009524E+03
#y = -1,2135365079E-04x3 + 4,8396812698E-02x2 - 6,3534816508E+00x + 2,7470335714E+02
#y = 4,1316721501E-04x3 - 1,6333499163E-01x2 + 2,1297012791E+01x - 9,1488577143E+02

"""
def rhoCH4(P, T, fluid):
    rho=PropsSI("D", "T", T, "P", P, fluid)
    return rho
def cpCH4(P,T, fluid):
    Cp=PropsSI("C", "P", P, "T", T, fluid)
    return Cp
def condCH4(P, T, fluid):
    conductivity=PropsSI("L", "T", T, "P", P, fluid)
    return conductivity
def viscCH4(P, T, fluid):
    viscosity=PropsSI("V", "T", T, "P", P, fluid)
    return viscosity
def sonCH4(P, T, fluid):
    celerite=PropsSI("A", "T", T, "P", P, fluid)
    return celerite

def entCH4(P,T,fluid):
    ent=PropsSI("S","T",T,"P",P,fluid)
    return ent
def pressureCH4(S,T,fluid):
    h1 = PropsSI('P', 'T',T, 'S',S,fluid)
    return h1

def DeltaT(P,T,fluid):
    temp=PropsSI("T", "P", P, "Q", 0, fluid)
    Delta_T=temp-T
    return Delta_T
"""
P = 3500000 #atm
T = 111 #K
Conductivity = condCH4(P, T)
Viscosity = viscCH4(P, T)
Density = rhoCH4(P, T)
Cp=cpCH4(P,T)
print(Conductivity, Viscosity, Density,Cp)
"""








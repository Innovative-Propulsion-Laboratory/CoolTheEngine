# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 18:44:22 2021

@author: julie
"""

def cas1(h,dx,lamb,Tf,inv):#validé
    a=1
    b=1
    c=0.5
    d=0.5
    x=-inv*(h*dx/lamb)-3
    plus=-(h*dx/lamb)*Tf
    return a,b,c,d,x,plus
def cas2(h,dx,lamb,Tf,inv):#validé
    a=0.5
    b=1
    c=0.5
    d=0
    x=-(h*dx/lamb)-2
    plus=-inv*(h*dx/lamb)*Tf#/((h*dx/Tf)+2)
    return a,b,c,d,x,plus
def cas3(h,dx,lamb,Tf,inv):#trivial validé
    a=0.25
    b=0.25
    c=0.25
    d=0.25
    x=-1
    plus=0
    return a,b,c,d,x,plus
def cas4(h,dx,lamb,Tf,inv):
    a=0.5
    b=0
    c=0
    d=0.5
    x=-inv*(1+(h*dx/(2*lamb)))
    plus=-(h*dx/(2*lamb))*Tf
    return a,b,c,d,x,plus
def cas5(h,dx,lamb,Tf,inv):
    h=0
    a=1
    b=0.5
    c=0
    d=0.5
    x=-(h*dx/lamb)-2
    plus=-inv*(h*dx/lamb)*Tf
    return a,b,c,d,x,plus

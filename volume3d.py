# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 17:37:58 2021

@author: julie
"""
from ProgressBar import ProgressBar
import numpy as np
from matplotlib import pyplot as plt


def carto3d(inv, x, y, mesure, col, title, number, limitation):
    if inv[0] == 1:
        x.reverse()
    if inv[1] == 1:
        y.reverse()
    if inv[2] == 1:
        mesure.reverse()
    b = 0
    Bar = ProgressBar(100, 30, "3D results computation          ")
    av = 100 / len(x)
    
    fig = plt.figure(figsize=(10, 17), dpi=500)
    ax = fig.add_subplot(111, projection='3d')
    Tcolor = []
    xu = []
    yu = []
    zu = []
    cu = []
    for i in range(0, len(x)):
        for k in range(0, number):
            theta = [k * (2 * np.pi / number)]
            for j in range(0, len(mesure[i]) - 1):
                the = theta[j] + (2 * np.pi / (number * len(mesure[i] * 2)))
                theta.append(the)
            # theta=np.linspace((2*np.pi/number)*k,(2*np.pi/number)+k*(2*np.pi/number),2*np.pi/(number*len(temp[0])))
            eh = 0
            for t in theta:
                yu.append(y[i] * np.cos(t))
                yu.append(y[i] * np.cos(-t))
                zu.append(y[i] * np.sin(t))
                zu.append(y[i] * np.sin(-t))
                xu.append(x[i])
                xu.append(x[i])
                cu.append(mesure[i][eh])
                cu.append(mesure[i][eh])
                eh += 1
        b += av
        Bar.update(b)
    print()
    Tcolor = [cu[i] for i in range(0, len(cu))]
    
    Bar = ProgressBar(100, 30, "3D results visualisation        ")
    Bar.update(10)
    p = ax.scatter(yu, zu, xu, c=Tcolor, marker='.', s=5, cmap=col)
    Bar.update(25)
    mis = min(xu)
    mas = max(xu)
    may = max(yu)
    ax.set_xlim(may + limitation, -may - limitation)
    ax.set_ylim(may + limitation, -may - limitation)
    ax.set_zlim(mas, mis)
    Bar.update(50)
    ax.view_init(15, 150)
    plt.title(title, fontsize=25)
    fig.colorbar(p, ax=ax, shrink=0.4, aspect=15)
    plt.show()
    Bar.update(75)
    if inv[0] == 1:
        x.reverse()
    if inv[1] == 1:
        y.reverse()   
    if inv[2] == 1:
        mesure.reverse()
    Bar.update(100)
    print()
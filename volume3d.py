# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 17:37:58 2021

@author: julie
"""
from ProgressBar import *
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def carto3d(inv, x, yprim, temp, col, title, number, limitation):
    if inv[0] == 1:
        x.reverse()
    if inv[1] == 1:
        yprim.reverse()
    if inv[2] == 1:
        temp.reverse()
    b = 0
    Bar = ProgressBar(100, 30, "Visualisation des résultats 3D  ")
    ay = 100 / (len(x))
    fig = plt.figure(figsize=(10, 17), dpi=500)
    ax = fig.add_subplot(111, projection='3d')
    Tcolor = []
    xu = []
    yu = []
    zu = []
    cu = []

    for i in range(0, len(x), 1):
        for k in range(0, number, 1):
            theta = [k * (2 * np.pi / number)]
            for j in range(0, len(temp[i]) - 1, 1):
                the = theta[j] + (2 * np.pi / (number * len(temp[i] * 2)))
                theta.append(the)
            # theta=np.linspace((2*np.pi/number)*k,(2*np.pi/number)+k*(2*np.pi/number),2*np.pi/(number*len(temp[0])))
            eh = 0
            for O in theta:
                yo = yprim[i] * np.cos(O)
                zo = yprim[i] * np.sin(O)
                yi = yprim[i] * np.cos(-O)
                zi = yprim[i] * np.sin(-O)
                xo = x[i]
                yu.append(yo)
                yu.append(yi)
                zu.append(zo)
                zu.append(zi)
                xu.append(xo)
                xu.append(xo)
                co = temp[i][eh]
                # print(i,eh)
                cu.append(co)
                cu.append(co)
                eh += 1
        b += ay
        Bar.update(b)

    mi = min(cu)
    ma = max(cu)
    for i in range(0, len(cu), 1):
        # T=[((cu[i]-min(cu))/max(cu))%1,(cu[i]-min(cu))/max(cu)%1,1-(cu[i]-min(cu))/max(cu)%1]
        # posey=(cu[i]-mi)/(ma-mi)

        # print(pos)
        T = cu[i]

        Tcolor.append(T)
        b += ay
        Bar.update(b)
    p = ax.scatter(yu, zu, xu, c=Tcolor, marker='.', s=5, cmap=col)
    mis = min(xu)
    mas = max(xu)
    miy = min(yu)
    may = max(yu)
    ax.set_xlim(may + limitation, -may - limitation)
    ax.set_ylim(may + limitation, -may - limitation)
    ax.set_zlim(mas, mis)
    ax.view_init(15, 150)
    plt.title(title, fontsize=25)
    fig.colorbar(p, ax=ax, shrink=0.4, aspect=15)
    plt.show()
    if inv[0] == 1:
        x.reverse()
    if inv[1] == 1:
        yprim.reverse()
    if inv[2] == 1:
        temp.reverse()
# col=plt.cm.hot
# view3d(x,yprim,temp,col)

# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 08:35:13 2021

@author: julie
"""

from ProgressBar import *
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def view3d(inv, x, yprim, temp, col, title, size2, limitation):
    if inv[0] == 1:
        x.reverse()
    if inv[1] == 1:
        yprim.reverse()
    if inv[2] == 1:
        temp.reverse()
    b = 0
    Bar = ProgressBar(100, 30, "Results visualization           ")
    ay = 100 / ((len(x) / 5) * 40 + len(x))
    fig = plt.figure(figsize=(10, size2), dpi=200)
    ax = fig.add_subplot(111, projection='3d')
    Tcolor = []
    xu = []
    yu = []
    zu = []
    cu = []

    for i in range(0, len(x), 2):
        theta = np.linspace(0, 2 * np.pi, 200)
        for O in theta:
            yo = yprim[i] * np.cos(O)
            zo = yprim[i] * np.sin(O)
            xo = x[i]
            yu.append(yo)
            zu.append(zo)
            xu.append(xo)
            co = temp[i]
            cu.append(co)
            b += ay
            Bar.update(b)
    print()
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
    p = ax.scatter(yu, zu, xu, c=Tcolor, marker='.', s=60, cmap=col)
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

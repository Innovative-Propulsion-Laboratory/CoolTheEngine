# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 08:35:13 2021

@author: julie
"""

from ProgressBar import ProgressBar
import numpy as np
from matplotlib import pyplot as plt


def view3d(inv, x, y, mesure, col, title, size2, limitation):
    if inv[0] == 1:
        x.reverse()
    if inv[1] == 1:
        y.reverse()
    if inv[2] == 1:
        mesure.reverse()

    fig = plt.figure(figsize=(10, size2), dpi=200)  # figsize=(float, float) : width, height in inches
    ax = fig.add_subplot(111, projection='3d')
    xu = []  # List of x position of points (vertical)
    yu = []  # List of y position of points
    zu = []  # List of z position of points
    cu = []  # List of mesure on each point
    theta = np.linspace(0, 2 * np.pi, 100)  # List of angles in order to make a whole circle

    for i in range(0, len(x), 5):  # step=3 to reduce computation time (this doesn't reduce quality much)
        for j in range(0, len(theta)):
            yu.append(y[i] * np.cos(theta[j]))
            zu.append(y[i] * np.sin(theta[j]))
            xu.append(x[i])
            cu.append(mesure[i])

    Tcolor = [cu[i] for i in range(0, len(cu))]
    p = ax.scatter(yu, zu, xu, c=Tcolor, marker='.', s=60, cmap=col)
    mis = min(xu)
    mas = max(xu)
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
        y.reverse()
    if inv[2] == 1:
        mesure.reverse()

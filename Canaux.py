# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 21:46:19 2020

@author: julien
"""
import numpy as np
import csv
import matplotlib.pyplot as plt


# x_coords_filename,y_coords_filename
# def canauxangl(plagex, plagey, nbc, lrg_col, ht, ht_c, ht_div, tore, debit_total, epaisseur_chemise, e_col, e_div, e_c):
#     figure_dpi = 150
#     crx = csv.reader(open(plagex, "r"))  # ouverture des x
#     cry = csv.reader(open(plagey, "r"))  # ouverture des y
#     x_value = []  # en m
#     y_value = []  # en m
#
#     for row in crx:  #
#         a = float(row[0]) / 1000  #
#         x_value.append(a)  # Récupération des données des deux fichiers
#     for row in cry:  #
#         a = float(row[0]) / 1000  #
#         y_value.append(a)  #
#     # print(x_value,y_value)
#
#     xcanauxre = []
#     ycanauxre = []
#     i = 0
#     while x_value[i] <= tore:
#         xcanauxre.append(x_value[i])
#         ycanauxre.append(y_value[i] + epaisseur_chemise)
#         i += 1
#     # print(xcanauxre,ycanauxre)
#     longcanal = len(xcanauxre)
#
#     degportion = 360 / nbc
#     poscol = ycanauxre.index(min(ycanauxre))
#     # print(poscol)
#     expcoef = (((2 * np.pi * ycanauxre[poscol]) / nbc) / lrg_col)
#     # print(expcoef)
#
#     larg_canalre = []
#     reste = []
#     for i in range(0, longcanal, 1):
#         larg = (ycanauxre[i] / ycanauxre[poscol]) * lrg_col
#         rest = ((ycanauxre[i] * 2 * np.pi) / nbc) - larg
#         larg_canalre.append(larg)
#         reste.append(rest)
#
#     n = 1
#     while ycanauxre[n] == ycanauxre[n - 1]:
#         n = n + 1
#     # print(n)
#     htre = []
#     for i in range(0, longcanal, 1):
#         if ycanauxre[i] == ycanauxre[0]:
#             htre.append(ht_c)
#         elif i <= poscol:
#             adcoef = (ht - ht_c) / (xcanauxre[poscol] - xcanauxre[n])
#             # print(adcoef)
#             hauteur = ht + adcoef * xcanauxre[i]
#             # print(hauteur)
#             htre.append(hauteur)
#         else:
#             augcoef = (ht_div - ht) / (xcanauxre[longcanal - 1] - xcanauxre[poscol])
#             hauteur = ht + augcoef * xcanauxre[i]
#             htre.append(hauteur)
#     # print(htre)
#     plt.figure(dpi=figure_dpi)
#     plt.plot(xcanauxre, ycanauxre, color='chocolate')
#     plt.title('trajet des canaux')
#     plt.show()
#     plt.figure(dpi=figure_dpi)
#     plt.plot(xcanauxre, larg_canalre, label='largeur', color='red')
#     plt.plot(xcanauxre, htre, label='hauteur', color='royalblue')
#     plt.plot(xcanauxre, reste, label='reste', color='chocolate')
#     plt.title('largeur, hauteur, reste des canaux')
#     plt.legend()
#     plt.show()
#
#     # calcul des aires
#     Areare = []
#     for i in range(0, longcanal, 1):
#         A = larg_canalre[i] * htre[i]
#         Areare.append(A)
#
#     # calcul indicatif des vitesses
#     vitessere = []
#     for i in range(0, longcanal, 1):
#         V = (debit_total / nbc) / Areare[i]
#         vitessere.append(V)
#     plt.figure(dpi=figure_dpi)
#     plt.plot(xcanauxre, Areare, color='chocolate')
#     plt.title('Aire des canaux')
#     plt.show()
#     plt.figure(dpi=figure_dpi)
#     plt.plot(xcanauxre, vitessere, color='chocolate')
#     plt.title('Vitesse représentative dans les canaux')
#     plt.show()
#     return xcanauxre, ycanauxre, larg_canalre, Areare, htre, reste


def canaux(profile_data, width_data, height_data, thickness_data, coefficients,
           tore_pos, debit_volumique_total_cool, nbc, plot_detail, write_in_csv, figure_dpi):
    """
    This function computes the caracteristics of channels on each point
    by interpolation between given values at injection plate (inj), end of cylindrical chamber (conv), 
    throat (col) and extremity of the nozzle (div).
    """

    x_value, y_value = profile_data
    lrg_inj, lrg_conv, lrg_col, lrg_tore = width_data
    ht_inj, ht_conv, ht_col, ht_tore = height_data
    e_conv, e_col, e_tore = thickness_data
    n1, n2, n3, n4, n5, n6 = coefficients

    xcanaux = []  # List of x where there are channels (before the manifold) (in m)
    y_coord_avec_canaux = []  # List of y where there are channels (before the manifold) (in m)
    # WARNING ! ycanaux will later in this function become y on the coolant side (in m)
    i = 0
    while i < len(x_value) and x_value[i] <= tore_pos:
        xcanaux.append(x_value[i])
        y_coord_avec_canaux.append(y_value[i])
        i += 1

    pos_conv = 0  # Index of the end of cylindrical chamber
    while y_coord_avec_canaux[pos_conv] == y_coord_avec_canaux[pos_conv + 1]:
        pos_conv += 1
    pos_col = y_coord_avec_canaux.index(min(y_coord_avec_canaux))  # Index of the throat
    y_coord_col = y_coord_avec_canaux[pos_col]  # y coordonate of the hot wall at the throat
    y_coord_inj = y_coord_avec_canaux[0]  # y coordonate of the hot wall at the injection plate
    y_coord_tore = y_coord_avec_canaux[-1]  # y coordonate of the hot wall at the manifold
    longc = len(xcanaux)  # Index of the end of the channels (at the manifold)

    wall_thickness = []  # Thickness of the chamber wall as a function of the engine axis (in m)
    acc = (e_conv - e_col) / (y_coord_inj - y_coord_col)
    for i in range(0, pos_col + 1):  # Chamber + convergent computation
        r = y_coord_avec_canaux[i]
        aug = (y_coord_col - r) / (y_coord_inj - y_coord_col)
        epp_x = ((-aug) ** n5) * (r - y_coord_col) * acc + e_col
        wall_thickness.append(epp_x)
    acc = (e_tore - e_col) / (y_coord_tore - y_coord_col)
    for i in range(pos_col + 1, longc):  # Divergent computation
        r = y_coord_avec_canaux[i]
        aug = (y_coord_tore - r) / (y_coord_tore - y_coord_col)
        epp_x = ((1 - aug) ** n6) * (r - y_coord_col) * acc + e_col
        wall_thickness.append(epp_x)

    angulaire = [0]
    ycanaux = [y_coord_inj + wall_thickness[0]]  # y of wall on coolant side (matched with engine geometry)
    for i in range(1, longc):
        vect = (xcanaux[i] - xcanaux[i - 1]) / ((((y_coord_avec_canaux[i] - y_coord_avec_canaux[i - 1]) ** 2) +
                                                  ((xcanaux[i] - xcanaux[i - 1]) ** 2)) ** 0.5)
        angulaire.append(np.rad2deg(np.arccos(vect)))
        """
        newep = ycanaux[i] + epaiss_chemise[i] / np.cos(np.deg2rad(angulaire[i]))
        ycanaux.append(newep)
        """
        ycanaux.append(y_coord_avec_canaux[i] + wall_thickness[i] / vect)

    if plot_detail >= 3:
        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, y_coord_avec_canaux, color='chocolate', label='y on hot gas side')
        plt.plot(xcanaux, ycanaux, color='blue', label='y on coolant side')
        plt.title('y coordinate of wall as a function of the engine axis')
        plt.legend(loc='lower left')
        plt.show()

    veritas = []
    for i in range(0, longc):
        verifepe = (((ycanaux[i] - y_value[i]) ** 2) - (
                np.sin(np.deg2rad(angulaire[i])) * (ycanaux[i] - y_value[i])) ** 2) ** 0.5
        veritas.append(verifepe)

    if plot_detail >= 3:
        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, veritas)
        plt.title('Channel thickness verification')
        plt.show()
        
        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, ycanaux, color='chocolate')
        plt.title('Channel travel as a function of the engine axis (y of the cold wall)')
        plt.show()

    debit_volumique_canal = debit_volumique_total_cool / nbc  # Volumic flow rate in a channel
    y_col = ycanaux[pos_col]  # y coordonate of the cold wall at the throat
    y_inj = ycanaux[0]  # y coordonate of the cold wall at the injection plate
    y_tore = ycanaux[-1]  # y coordonate of the cold wall at the manifold
    
    larg_ailette = []  # Width of a rib as a function of the engine axis (in m)
    larg_canal = []  # Width of a channel as a function of the engine axis (in m)
    pente = (lrg_conv - lrg_inj) / (xcanaux[pos_conv] - xcanaux[0])
    for i in range(0, pos_conv + 1):  # Chamber computation
        r = ycanaux[i]
        lrg_x = pente * (xcanaux[i] - xcanaux[0]) + lrg_inj
        lrg_aill = (r * 2 * np.pi / nbc) - lrg_x
        larg_ailette.append(lrg_aill)
        larg_canal.append(lrg_x)
    acc = (lrg_conv - lrg_col) / (y_inj - y_col)
    for i in range(pos_conv + 1, pos_col + 1):  # Convergent computation
        r = ycanaux[i]
        aug = (y_col - r) / (y_inj - y_col)
        lrg_x = ((-aug) ** n1) * (r - y_col) * acc + lrg_col
        lrg_aill = (r * 2 * np.pi / nbc) - lrg_x
        larg_ailette.append(lrg_aill)
        larg_canal.append(lrg_x)
    acc = (lrg_tore - lrg_col) / (y_tore - y_col)
    for i in range(pos_col + 1, longc):  # Divergent computation
        r = ycanaux[i]
        aug = (y_tore - r) / (y_tore - y_col)
        lrg_x = ((1 - aug) ** n2) * (r - y_col) * acc + lrg_col
        lrg_aill = (r * 2 * np.pi / nbc) - lrg_x
        larg_ailette.append(lrg_aill)
        larg_canal.append(lrg_x)

    ht_canal = []  # Height of a channel as a function of the engine axis (in m)
    pente = (ht_conv - ht_inj) / (xcanaux[pos_conv] - xcanaux[0])
    for i in range(0, pos_conv + 1):  # Chamber computation
        htr_x = pente * (xcanaux[i] - xcanaux[0]) + ht_inj
        ht_canal.append(htr_x)
    acc = (ht_conv - ht_col) / (y_inj - y_col)
    for i in range(pos_conv + 1, pos_col):  # Convergent computation
        r = ycanaux[i]
        aug = (y_col - r) / (y_inj - y_col)
        htr_x = ((-aug) ** n3) * (r - y_col) * acc + ht_col
        ht_canal.append(htr_x)
    acc = (ht_tore - ht_col) / (y_tore - y_col)
    for i in range(pos_col, longc):  # Divergent computation
        r = ycanaux[i]
        aug = (y_tore - r) / (y_tore - y_col)
        htr_x = ((1 - aug) ** n4) * (r - y_col) * acc + ht_col
        ht_canal.append(htr_x)
    """
    for zz in larg_canal:
     #A supprimer
        x=xcanaux[larg_canal.index(zz)]
        z2=zz+0.5*zz*sin(x*1000/4)
        larg_canal[larg_canal.index(zz)]=z2 
    """
    area_channel = []  # Area of a channel as a function of the engine axis (m²)
    vitesse_coolant = []  # Velocity of coolant in a channel as a function of the engine axis (m/s)
    for i in range(0, longc):
        aire = larg_canal[i] * ht_canal[i]
        area_channel.append(aire)
        v = debit_volumique_canal / aire
        vitesse_coolant.append(v)

    if write_in_csv:
        "Writing the results of the study in a CSV file"
        file_name = "output/channel_macro_catia.csv"
        file = open(file_name, "w", newline="")
        writer = csv.writer(file)
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow((1000 * xcanaux[i], 1000 * (ycanaux[i]), 1000 * (larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow((1000 * xcanaux[i], 1000 * (ycanaux[i]), 1000 * (-larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow((1000 * xcanaux[i], 1000 * (ycanaux[i] + ht_canal[i]), 1000 * (larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow((1000 * xcanaux[i], 1000 * (ycanaux[i] + ht_canal[i]), 1000 * (- larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["End"])
        file.close()

    if plot_detail >= 3:
        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, larg_ailette, label='Rib width', color='chocolate')
        plt.plot(xcanaux, larg_canal, label='Channel width', color='green')
        plt.plot(xcanaux, ht_canal, label='Channel height', color='blue')
        plt.title('Width of channels and ribs')
        plt.legend(loc='upper left')
        plt.show()

        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, wall_thickness, color='chocolate')
        plt.title('Thickness of chamber wall as a function of the engine axis')
        plt.show()

        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, ht_canal, color='chocolate')
        plt.title('Channel height as a function of the engine axis')
        plt.show()

    if plot_detail >= 1:
        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, area_channel, color='chocolate')
        plt.title('Channel cross-sectionnal area as a function of the engine axis')
        plt.show()

    return xcanaux, ycanaux, larg_canal, larg_ailette, ht_canal, wall_thickness, area_channel, longc, y_coord_avec_canaux

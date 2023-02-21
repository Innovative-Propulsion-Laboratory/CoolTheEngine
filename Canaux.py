# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 21:46:19 2020

@author: julien
"""
import numpy as np
import csv
import matplotlib.pyplot as plt

"""
# x_coords_filename,y_coords_filename
def canauxangl(plagex, plagey, nbc, lrg_col, ht, ht_c, ht_div, tore, debit_total, epaisseur_chemise, e_col, e_div, e_c):
    crx = csv.reader(open(plagex, "r"))  # ouverture des x
    cry = csv.reader(open(plagey, "r"))  # ouverture des y
    x_value = []  # en m
    y_value = []  # en m

    for row in crx:  #
        a = float(row[0]) / 1000  #
        x_value.append(a)  # Récupération des données des deux fichiers
    for row in cry:  #
        a = float(row[0]) / 1000  #
        y_value.append(a)  #
    # print(x_value,y_value)

    xcanauxre = []
    ycanauxre = []
    i = 0
    while x_value[i] <= tore:
        xcanauxre.append(x_value[i])
        ycanauxre.append(y_value[i] + epaisseur_chemise)
        i += 1
    # print(xcanauxre,ycanauxre)
    longcanal = len(xcanauxre)

    degportion = 360 / nbc
    poscol = ycanauxre.index(min(ycanauxre))
    # print(poscol)
    expcoef = (((2 * np.pi * ycanauxre[poscol]) / nbc) / lrg_col)
    # print(expcoef)

    larg_canalre = []
    reste = []
    for i in range(0, longcanal, 1):
        larg = (ycanauxre[i] / ycanauxre[poscol]) * lrg_col
        rest = ((ycanauxre[i] * 2 * np.pi) / nbc) - larg
        larg_canalre.append(larg)
        reste.append(rest)

    n = 1
    while ycanauxre[n] == ycanauxre[n - 1]:
        n = n + 1
    # print(n)
    htre = []
    for i in range(0, longcanal, 1):
        if ycanauxre[i] == ycanauxre[0]:
            htre.append(ht_c)
        elif i <= poscol:
            adcoef = (ht - ht_c) / (xcanauxre[poscol] - xcanauxre[n])
            # print(adcoef)
            hauteur = ht + adcoef * xcanauxre[i]
            # print(hauteur)
            htre.append(hauteur)
        else:
            augcoef = (ht_div - ht) / (xcanauxre[longcanal - 1] - xcanauxre[poscol])
            hauteur = ht + augcoef * xcanauxre[i]
            htre.append(hauteur)
    # print(htre)
    plt.figure(dpi=200)
    plt.plot(xcanauxre, ycanauxre, color='chocolate')
    plt.title('trajet des canaux')
    plt.show()
    plt.figure(dpi=200)
    plt.plot(xcanauxre, larg_canalre, label='largeur', color='red')
    plt.plot(xcanauxre, htre, label='hauteur', color='royalblue')
    plt.plot(xcanauxre, reste, label='reste', color='chocolate')
    plt.title('largeur, hauteur, reste des canaux')
    plt.legend()
    plt.show()

    # calcul des aires
    Areare = []
    for i in range(0, longcanal, 1):
        A = larg_canalre[i] * htre[i]
        Areare.append(A)

    # calcul indicatif des vitesses
    vitessere = []
    for i in range(0, longcanal, 1):
        V = (debit_total / nbc) / Areare[i]
        vitessere.append(V)
    plt.figure(dpi=200)
    plt.plot(xcanauxre, Areare, color='chocolate')
    plt.title('Aire des canaux')
    plt.show()
    plt.figure(dpi=200)
    plt.plot(xcanauxre, vitessere, color='chocolate')
    plt.title('Vitesse représentative dans les canaux')
    plt.show()
    return xcanauxre, ycanauxre, larg_canalre, Areare, htre, reste
"""


def canaux(x_value, y_value, nbc, lrg_inj, lrg_conv, lrg_col, lrg_tore, ht_inj, ht_conv, ht_col, ht_tore,
           e_conv, e_col, e_tore, tore, debit_total, n1, n2, n3, n4, n5, n6):
    """
    This function compute the caracteristics of channels on each point 
    by interpolation between given values at injection plate (inj), end of cylindrical chamber (conv), 
    throat (col) and extremity of the nozzle (div).
    """
    xcanauxre = []  # List of x where there are channels (before the manifold) (in m)
    ycanauxre = []  # List of y where there are channels (before the manifold) (in m)
    i = 0
    while i < len(x_value) and x_value[i] <= tore:
        xcanauxre.append(x_value[i])
        ycanauxre.append(y_value[i])
        i += 1

    pos_conv = 0  # Index of the end of cylindrical chamber
    while ycanauxre[pos_conv] == ycanauxre[pos_conv + 1]:
        pos_conv += 1
    pos_col = ycanauxre.index(min(ycanauxre))  # Index of the throat
    y_col = ycanauxre[pos_col]  # y coordonate of the hot wall at the throat
    y_inj = ycanauxre[0]  # y coordonate of the hot wall at the injection plate
    y_tore = ycanauxre[-1]  # y coordonate of the hot wall at the manifold
    longc = len(xcanauxre)  # Index of the end of the channels (at the manifold)

    epaiss_chemise = []  # Thickness of the chamber wall as a function of the engine axis (in m)
    acc = (e_conv - e_col) / (y_inj - y_col)
    for i in range(0, pos_col + 1):  # Chamber + convergent computation
        r = ycanauxre[i]
        aug = ((y_col - r) / (y_inj - y_col))
        epp_x = ((-aug) ** n5) * (r - y_col) * acc + e_col
        epaiss_chemise.append(epp_x)
    acc = (e_tore - e_col) / (y_tore - y_col)
    for i in range(pos_col + 1, longc):  # Divergent computation
        r = ycanauxre[i]
        aug = ((y_tore - r) / (y_tore - y_col))
        epp_x = ((1 - aug) ** n6) * (r - y_col) * acc + e_col
        epaiss_chemise.append(epp_x)

    angulaire = [0]
    newepaisseur = [y_inj + epaiss_chemise[0]]  # Corrected thickness (to match with the geometry of the engine)
    for i in range(1, longc):
        vect2 = (xcanauxre[i] - xcanauxre[i - 1]) / ((((ycanauxre[i] - ycanauxre[i - 1]) ** 2) +
                                                      ((xcanauxre[i] - xcanauxre[i - 1]) ** 2)) ** 0.5)
        angulaire.append(np.rad2deg(np.arccos(vect2)))
        newep = ycanauxre[i] + epaiss_chemise[i] / np.cos(np.deg2rad(angulaire[i]))
        newepaisseur.append(newep)

    plt.plot(xcanauxre, ycanauxre, color='chocolate', label='old')  # (configuration qui
    plt.plot(xcanauxre, newepaisseur, color='blue', label='New')
    plt.title('Verif')  # et au Viserion)
    plt.legend()
    plt.show()

    ycanauxre = [newepaisseur[i] for i in range(0, len(newepaisseur))]
    veritas = []
    for i in range(0, longc):
        verifepe = (((ycanauxre[i] - y_value[i]) ** 2) - (
                np.sin(np.deg2rad(angulaire[i])) * (ycanauxre[i] - y_value[i])) ** 2) ** 0.5
        veritas.append(verifepe)
    """
    plt.plot(xcanauxre, veritas)  # (configuration qui
    plt.title('vérification epaisseur canaux')  # et au Viserion)
    plt.show()

    plt.plot(xcanauxre, ycanauxre, color='chocolate')  # (configuration qui
    plt.title('trajet des canaux')  # et au Viserion)
    plt.show()
    """
    debitcanal = debit_total / nbc  # Volumic debit of a channel
    y_col = ycanauxre[pos_col]  # y coordonate of the cold wall at the throat
    y_inj = ycanauxre[0]  # y coordonate of the cold wall at the injection plate
    y_tore = ycanauxre[-1]  # y coordonate of the cold wall at the manifold
    larg_ailette = []  # Width of a rib as a function of the engine axis (in m)
    larg_canalre = []  # Width of a channel as a function of the engine axis (in m)
    pente = (lrg_conv - lrg_inj) / (xcanauxre[pos_conv] - xcanauxre[0])
    for i in range(0, pos_conv + 1):  # Chamber computation
        r = ycanauxre[i]
        lrg_x = pente * (xcanauxre[i] - xcanauxre[0]) + lrg_inj
        lrg_aill = (r * 2 * np.pi / nbc) - lrg_x
        larg_ailette.append(lrg_aill)
        larg_canalre.append(lrg_x)
    acc = (lrg_conv - lrg_col) / (y_inj - y_col)
    for i in range(pos_conv + 1, pos_col + 1):  # Convergent computation
        r = ycanauxre[i]
        aug = ((y_col - r) / (y_inj - y_col))
        lrg_x = ((-aug) ** n1) * (r - y_col) * acc + lrg_col
        lrg_aill = (r * 2 * np.pi / nbc) - lrg_x
        larg_ailette.append(lrg_aill)
        larg_canalre.append(lrg_x)
    acc = (lrg_tore - lrg_col) / (y_tore - y_col)
    for i in range(pos_col + 1, longc):  # Divergent computation
        r = ycanauxre[i]
        aug = ((y_tore - r) / (y_tore - y_col))
        lrg_x = ((1 - aug) ** n2) * (r - y_col) * acc + lrg_col
        lrg_aill = (r * 2 * np.pi / nbc) - lrg_x
        larg_ailette.append(lrg_aill)
        larg_canalre.append(lrg_x)

    htre = []  # Height of a channel as a function of the engine axis (in m)
    pente = (ht_conv - ht_inj) / (xcanauxre[pos_conv] - xcanauxre[0])
    for i in range(0, pos_conv + 1):  # Chamber computation
        htr_x = pente * (xcanauxre[i] - xcanauxre[0]) + ht_inj
        htre.append(htr_x)
    acc = (ht_conv - ht_col) / (y_inj - y_col)
    for i in range(pos_conv + 1, pos_col):  # Convergent computation
        r = ycanauxre[i]
        aug = ((y_col - r) / (y_inj - y_col))
        htr_x = ((-aug) ** n3) * (r - y_col) * acc + ht_col
        htre.append(htr_x)
    acc = (ht_tore - ht_col) / (y_tore - y_col)
    for i in range(pos_col, longc):  # Divergent computation
        r = ycanauxre[i]
        aug = ((y_tore - r) / (y_tore - y_col))
        htr_x = ((1 - aug) ** n4) * (r - y_col) * acc + ht_col
        htre.append(htr_x)
    """
    for zz in larg_canalre:
     #A supprimer
        x=xcanauxre[larg_canalre.index(zz)]
        z2=zz+0.5*zz*sin(x*1000/4)
        larg_canalre[larg_canalre.index(zz)]=z2 
    """
    Areare = []  # Area of a channel as a function of the engine axis (without taking changes of density in count) (in m²)
    vitessere = []  # Velocity of coolant in a channel as a function of the engine axis (in m/s)
    for i in range(0, len(larg_canalre)):
        aire = larg_canalre[i] * htre[i]
        Areare.append(aire)
        v = debitcanal / aire
        vitessere.append(v)

    "Writing the results of the study in a CSV file"
    file_name = "channel_macro_catia.csv"
    file = open(file_name, "w", newline="")
    writer = csv.writer(file)
    writer.writerow(["StartCurve"])
    for i in range(0, longc, 3):
        writer.writerow((1000 * xcanauxre[i], 1000 * (ycanauxre[i]), 1000 * (larg_canalre[i] / 2)))
    writer.writerow(["EndCurve"])
    writer.writerow(["StartCurve"])
    for i in range(0, longc, 3):
        writer.writerow((1000 * xcanauxre[i], 1000 * (ycanauxre[i]), 1000 * (-larg_canalre[i] / 2)))
    writer.writerow(["EndCurve"])
    writer.writerow(["StartCurve"])
    for i in range(0, longc, 3):
        writer.writerow((1000 * xcanauxre[i], 1000 * (ycanauxre[i] + htre[i]), 1000 * (larg_canalre[i] / 2)))
    writer.writerow(["EndCurve"])
    writer.writerow(["StartCurve"])
    for i in range(0, longc, 3):
        writer.writerow((1000 * xcanauxre[i], 1000 * (ycanauxre[i] + htre[i]), 1000 * (- larg_canalre[i] / 2)))
    writer.writerow(["EndCurve"])
    writer.writerow(["End"])
    file.close()

    plt.plot(xcanauxre, larg_ailette, label='Rib width', color='chocolate')
    plt.plot(xcanauxre, larg_canalre, label='Channel width', color='green')
    plt.plot(xcanauxre, htre, label='Channel height', color='blue')
    plt.title('Width of channels and ribs')
    plt.legend()
    plt.show()

    plt.plot(xcanauxre, vitessere, color='chocolate')
    plt.title('Velocity of coolant in channels (in m/s) as a function of the engine axis')
    plt.show()

    plt.plot(xcanauxre, epaiss_chemise, color='chocolate')
    plt.title('Thickness of chamber wall as a function of the engine axis')
    plt.show()

    plt.plot(xcanauxre, htre, color='chocolate')
    plt.title('Channel height as a function of the engine axis')
    plt.show()

    plt.plot(xcanauxre, Areare, color='chocolate')
    plt.title('Channel cross-sectionnal area as a function of the engine axis')
    plt.show()

    return xcanauxre, ycanauxre, larg_canalre, larg_ailette, htre, epaiss_chemise, Areare, longc

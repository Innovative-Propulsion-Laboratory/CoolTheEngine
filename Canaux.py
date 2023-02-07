# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 21:46:19 2020

@author: julien
"""
import numpy as np
import csv
import matplotlib.pyplot as plt


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


def canaux(plagex, plagey, nbc, lrg, lrg_c, lrg_div, ht, ht_c, ht_div, tore, debit_total, n1, n2, n3, n4, e_col, e_div,
           e_c, n5, n6, lrg_c2, ht_c2):
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

    i = 0
    xcanauxre = []  # allé
    ycanauxre = []
    while x_value[i] <= tore:
        xcanauxre.append(x_value[i])
        ycanauxre.append(y_value[i])
        i += 1
    section = 0
    while y_value[section] == y_value[section + 1]:
        section += 1
    """
    angles=[]
    newxhtre=[]
    newyhtre=[]
    for i in range (0,len(xcanauxre),1):
        if i==0:
            angle=0
            angles.append(angle)
        elif i==(len(xcanauxre)-1):
            angle=angles[i-1]
            angles.append(angle)
        else:
            vect1=(xcanauxre[i]-xcanauxre[i-1])/((((ycanauxre[i]-ycanauxre[i-1])**2)+((xcanauxre[i]-xcanauxre[
            i-1])**2))**0.5)
            vect2=(xcanauxre[i+1]-xcanauxre[i])/((((ycanauxre[i+1]-ycanauxre[i])**2)+((xcanauxre[i+1]-xcanauxre[
            i])**2))**0.5)
            angle1=degrees(acos(vect1))
            angle2=degrees(acos(vect2))
            angle=angle2
            angles.append(angle)
        
        newx=xcanauxre[i]+htre[i]*sin(radians(angles[i]))
        newy=ycanauxre[i]+htre[i]*cos(radians(angles[i]))
        newxhtre.append(newx)
        newyhtre.append(newy)
    """
    col = ycanauxre.index(min(ycanauxre))  # Définition position du col (index)
    angulaire = []
    newepaisseur = []
    epaiss_chemise = []
    for i in range(0, len(xcanauxre), 1):
        x = xcanauxre[i]
        pos = xcanauxre.index(x)  # Indexation de la
        r = ycanauxre[pos]  # position considérée
        p = np.pi * 2 * r / nbc  #

        if (ycanauxre.index(ycanauxre[pos])) > col:  # Calcul largeur transition
            acc = (e_div - e_col) / (ycanauxre[-1] - ycanauxre[col])
            aug = ((ycanauxre[-1] - ycanauxre[pos]) / (ycanauxre[-1] - ycanauxre[col]))
            epp_x = ((1 - aug) ** n5) * (r - ycanauxre[col]) * acc + e_col
        else:
            acc = (e_c - e_col) / (ycanauxre[0] - ycanauxre[col])  # Calcul largeur canaux
            aug = ((ycanauxre[col] - ycanauxre[pos]) / (ycanauxre[0] - ycanauxre[col]))
            epp_x = ((-aug) ** n6) * (r - ycanauxre[col]) * acc + e_col

        epaiss_chemise.append(epp_x)

    for i in range(0, len(xcanauxre), 1):
        if i == 0:
            angle = 0
            angulaire.append(angle)
        elif i == (len(xcanauxre) - 1):
            vect2 = (xcanauxre[i] - xcanauxre[i - 1]) / (
                    (((ycanauxre[i] - ycanauxre[i - 1]) ** 2) + ((xcanauxre[i] - xcanauxre[i - 1]) ** 2)) ** 0.5)
            angle2 = np.rad2deg(np.arccos(vect2))
            angulaire.append(angle)
        else:
            vect2 = (xcanauxre[i] - xcanauxre[i - 1]) / (
                    (((ycanauxre[i] - ycanauxre[i - 1]) ** 2) + ((xcanauxre[i] - xcanauxre[i - 1]) ** 2)) ** 0.5)
            angle2 = np.rad2deg(np.arccos(vect2))
            angle = angle2
            angulaire.append(angle)

        newep = ycanauxre[i] + epaiss_chemise[i] / np.cos(np.deg2rad(angulaire[i]))
        newepaisseur.append(newep)
    ycanauxre = []
    for i in range(0, len(newepaisseur), 1):
        ycanauxre.append(newepaisseur[i])
    veritas = []
    for i in range(0, len(xcanauxre), 1):
        verifepe = (((ycanauxre[i] - y_value[i]) ** 2) - (
                np.sin(np.deg2rad(angulaire[i])) * (ycanauxre[i] - y_value[i])) ** 2) ** 0.5
        veritas.append(verifepe)

    # plt.plot(xcanauxre, veritas)  # (configuration qui
    # plt.title('vérification epaisseur canaux')  # et au Viserion)
    # plt.show()

    # plt.plot(xcanauxre, ycanauxre, color='chocolate')  # (configuration qui
    # plt.title('trajet des canaux')  # et au Viserion)
    # plt.show()
    ###################################################
    debitcanal = debit_total / nbc  # calcul débit volumique par canaux

    larg_ailette = []  #
    larg_canalre = []  #

    for x in xcanauxre:  #
        pos = xcanauxre.index(x)  # Indexation de la
        r = ycanauxre[pos]  # position considérée
        p = np.pi * 2 * r / nbc  #
        # print(xcanauxre[section])
        # print(xcanauxre[0])
        pente = (lrg_c - lrg_c2) / (xcanauxre[section] - xcanauxre[0])
        if (ycanauxre.index(ycanauxre[pos])) > col:  # Calcul largeur transition
            acc = (lrg_div - lrg) / (ycanauxre[-1] - ycanauxre[col])
            aug = ((ycanauxre[-1] - ycanauxre[pos]) / (ycanauxre[-1] - ycanauxre[col]))
            lrg_x = ((1 - aug) ** n2) * (r - ycanauxre[col]) * acc + lrg
            lrg_aill = (r * 2 * np.pi / nbc) - lrg_x
        else:
            if pos <= section:
                lrg_x = pente * (x - xcanauxre[0]) + lrg_c2
            else:
                acc = (lrg_c - lrg) / (ycanauxre[0] - ycanauxre[col])  # Calcul largeur canaux
                aug = ((ycanauxre[col] - ycanauxre[pos]) / (ycanauxre[0] - ycanauxre[col]))
                lrg_x = ((-aug) ** n1) * (r - ycanauxre[col]) * acc + lrg
            lrg_aill = (r * 2 * np.pi / nbc) - lrg_x
        larg_ailette.append(lrg_aill)
        larg_canalre.append(lrg_x)

    htre = []
    for x in xcanauxre:  #
        pos = xcanauxre.index(x)  # Indexation de la
        r = ycanauxre[pos]  # position considérée
        p = np.pi * 2 * r / nbc  #
        pente = (ht_c - ht_c2) / (xcanauxre[section] - xcanauxre[0])
        if (xcanauxre[pos]) >= 0:  # Calcul largeur transition
            acc = (ht_div - ht) / (ycanauxre[-1] - ycanauxre[col])
            aug = ((ycanauxre[-1] - ycanauxre[pos]) / (ycanauxre[-1] - ycanauxre[col]))
            htr_x = ((1 - aug) ** n4) * (r - ycanauxre[col]) * acc + ht

        else:
            if pos <= section:
                htr_x = pente * (x - xcanauxre[0]) + ht_c2
            else:
                acc = (ht_c - ht) / (ycanauxre[0] - ycanauxre[col])  # Calcul largeur canaux
                aug = ((ycanauxre[col] - ycanauxre[pos]) / (ycanauxre[0] - ycanauxre[col]))
                htr_x = ((-aug) ** n3) * (r - ycanauxre[col]) * acc + ht  # sur le reste du moteur

        htre.append(htr_x)
    """
    for zz in larg_canalre:
     #A supprimer
        x=xcanauxre[larg_canalre.index(zz)]
        z2=zz+0.5*zz*sin(x*1000/4)
        larg_canalre[larg_canalre.index(zz)]=z2 
    """

    Areare = []  # (sans prendre en compte
    vitessere = []  #
    rang = 0  # changement de rho)
    for y in larg_canalre:  #
        aire = y * htre[rang]
        Areare.append(aire)  #
        v = debitcanal / aire  #
        vitessere.append(v)
        rang += 1

        #

    plt.plot(xcanauxre, larg_canalre, label='Canal width', color='green')  # Création et affichage
    plt.plot(xcanauxre, larg_ailette, label='Fin width', color='chocolate')
    plt.plot(xcanauxre, htre, label='Height', color='blue')  #
    plt.title('Width of channels and fins (ailettes in french)')  # bleu=retour
    plt.legend()
    plt.show()  # orange=allé
    plt.plot(xcanauxre, vitessere, color='chocolate')  # Affichage des courbes
    plt.title('Velocity of coolant in channels (in m/s) as a function of engine axis')  # bleu=retour
    plt.show()

    plt.plot(xcanauxre, epaiss_chemise, color='chocolate')  # Affichage des courbes
    plt.title('Thickness of chamber wall as a function of engine axis')  # bleu=retour
    plt.show()

    plt.plot(xcanauxre, htre, color='chocolate')
    plt.title('Channel height as a function of engine axis')  # bleu=retour
    plt.show()  # orange=allé

    plt.plot(xcanauxre, Areare, color='chocolate')
    plt.title('Channel cross-sectionnal area as a function of engine axis')  # bleu=retour
    plt.show()

    return xcanauxre, ycanauxre, larg_canalre, Areare, htre, larg_ailette, epaiss_chemise

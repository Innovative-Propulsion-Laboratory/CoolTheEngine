# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 21:46:19 2020

@author: julien
"""
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
import cte_tools as t


def canaux(profile_data, width_data, height_data, thickness_data, coefficients,
           tore_pos, debit_volumique_total_cool, nbc, plot_detail, write_in_csv, figure_dpi, plot_dir=None):
    """
    This function computes the caracteristics of channels on each point
    by interpolation between given values at injection plate (inj), end of cylindrical chamber (conv), 
    throat (col) and extremity of the nozzle (div).
    """
    if plot_detail >= 1 and plot_dir is not None:
        plot_dir = os.path.join(plot_dir, "Channels")
        if not os.path.exists(plot_dir) and not os.path.isdir(plot_dir):
            os.mkdir(plot_dir)

    x_value, y_value = profile_data
    lrg_inj, lrg_conv, lrg_col, lrg_tore = width_data
    ht_inj, ht_conv, ht_col, ht_tore = height_data
    e_conv, e_col, e_tore = thickness_data
    n1, n2, n3, n4, n5, n6 = coefficients

    # List of x where there are channels (before the manifold) (in m)
    xcanaux = []
    # List of y where there are channels (before the manifold) (in m)
    y_coord_avec_canaux = []
    i = 0
    while i < len(x_value) and x_value[i] <= tore_pos:
        xcanaux.append(x_value[i])
        y_coord_avec_canaux.append(y_value[i])
        i += 1
    xcanaux_mm = [xc * 1000 for xc in xcanaux]

    pos_conv = 0  # Index of the end of cylindrical chamber
    while y_coord_avec_canaux[pos_conv] == y_coord_avec_canaux[pos_conv + 1]:
        pos_conv += 1
    pos_col = y_coord_avec_canaux.index(
        min(y_coord_avec_canaux))  # Index of the throat
    # y coordonate of the hot wall at the throat
    y_coord_col = y_coord_avec_canaux[pos_col]
    # y coordonate of the hot wall at the injection plate
    y_coord_inj = y_coord_avec_canaux[0]
    # y coordonate of the hot wall at the manifold
    y_coord_tore = y_coord_avec_canaux[-1]
    longc = len(xcanaux)  # Index of the end of the channels (at the manifold)

    # Thickness of the chamber wall as a function of the engine axis (in m)
    wall_thickness = []
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
    # y of wall on coolant side (matched with engine geometry)
    ycanaux = [y_coord_inj + wall_thickness[0]]
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
        fig = t.n_plots(x=xcanaux, y_list=[y_coord_avec_canaux, ycanaux],
                        xlabel=r'Position x along the engine [m]',
                        ylabel=r'Radius [m]',
                        y_label_list=["Gas-side", "Coolant-side"],
                        colors_list=["red", "blue"],
                        title="y coordinate of the wall",
                        show=False, xmin=-0.190, xmax=0.035, ymin=0)
        if plot_dir is not None:
            fig.savefig(f"{plot_dir}/Wall y.png")
            plt.close(fig)
        else:
            plt.show()

    veritas = []
    for i in range(0, longc):
        verifepe = (((ycanaux[i] - y_value[i]) ** 2) - (
                np.sin(np.deg2rad(angulaire[i])) * (ycanaux[i] - y_value[i])) ** 2) ** 0.5
        veritas.append(verifepe)

    if plot_detail >= 3:
        fig = t.one_plot(x=xcanaux, y=veritas,
                         xlabel=r'Position x along the engine [m]',
                         ylabel=r'Distance [m]',
                         title="Channel thickness verification",
                         show=False, xmin=-0.190, xmax=0.035, sci_notation=True)
        if plot_dir is not None:
            fig.savefig(f"{plot_dir}/Channel thickness verif.png")
            plt.close(fig)
        else:
            plt.show()

        fig = t.one_plot(x=xcanaux, y=ycanaux,
                         xlabel=r'Position x along the engine [m]',
                         ylabel=r'Radius [m]',
                         title="Y of the cold wall",
                         show=False, xmin=-0.190, xmax=0.035,
                         ymin=0, sci_notation=True)
        if plot_dir is not None:
            fig.savefig(f"{plot_dir}/Channel travel.png")
            plt.close(fig)
        else:
            plt.show()

    debit_volumique_canal = debit_volumique_total_cool / \
                            nbc  # Volumic flow rate in a channel
    y_col = ycanaux[pos_col]  # y coordonate of the cold wall at the throat
    y_inj = ycanaux[0]  # y coordonate of the cold wall at the injection plate
    y_tore = ycanaux[-1]  # y coordonate of the cold wall at the manifold

    larg_ailette = []  # Width of a rib as a function of the engine axis (in m)
    # Width of a channel as a function of the engine axis (in m)
    larg_canal = []
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

    # Height of a channel as a function of the engine axis (in m)
    ht_canal = []
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

    area_channel = []  # Area of a channel as a function of the engine axis (m²)
    # Velocity of coolant in a channel as a function of the engine axis (m/s)
    vitesse_coolant = []
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
            writer.writerow(
                (1000 * xcanaux[i], 1000 * (ycanaux[i]), 1000 * (larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow(
                (1000 * xcanaux[i], 1000 * (ycanaux[i]), 1000 * (-larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow(
                (1000 * xcanaux[i], 1000 * (ycanaux[i] + ht_canal[i]), 1000 * (larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow(
                (1000 * xcanaux[i], 1000 * (ycanaux[i] + ht_canal[i]), 1000 * (- larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["End"])
        file.close()

    if plot_detail >= 3:
        fig = t.n_plots(x=xcanaux,
                        y_list=[larg_ailette,
                                larg_canal,
                                [1e-3] * len(larg_canal)],
                        xlabel=r'Position x along the engine [m]',
                        ylabel=r'Width [m]',
                        y_label_list=["Rib width", "Channel width", "1 mm"],
                        colors_list=["chocolate", "green", "red"],
                        title="Width of the channels and ribs",
                        show=False, xmin=-0.190, xmax=0.035,
                        ymin=0, sci_notation=True)
        if plot_dir is not None:
            fig.savefig(f"{plot_dir}/Width of channels and ribs.png")
            plt.close(fig)
        else:
            plt.show()

        fig = t.one_plot(x=xcanaux, y=wall_thickness,
                         xlabel=r'Position x along the engine [m]',
                         ylabel=r'Thickness [m]',
                         title="Thickness of the chamber wall",
                         show=False, xmin=-0.190, xmax=0.035, sci_notation=True)
        if plot_dir is not None:
            fig.savefig(f"{plot_dir}/Thickness of chamber wall.png")
            plt.close(fig)
        else:
            plt.show()

        fig = t.one_plot(x=xcanaux, y=ht_canal,
                         xlabel=r'Position x along the engine [m]',
                         ylabel=r'Height [m]',
                         title="Channel height",
                         show=False, xmin=-0.190, xmax=0.035,
                         ymin=0, sci_notation=True)
        if plot_dir is not None:
            fig.savefig(f"{plot_dir}/Channel height.png")
            plt.close(fig)
        else:
            plt.show()

    if plot_detail >= 1:
        fig = t.one_plot(x=xcanaux, y=area_channel,
                         xlabel=r'Position x along the engine [m]',
                         ylabel=r'Area [m²]',
                         title="Channel cross-sectionnal area",
                         show=False, xmin=-0.190, xmax=0.035,
                         ymin=0, sci_notation=True)
        if plot_dir is not None:
            fig.savefig(f"{plot_dir}/Channel area.png")
            plt.close(fig)
        else:
            plt.show()

    return xcanaux, ycanaux, larg_canal, larg_ailette, ht_canal, wall_thickness, area_channel, longc, y_coord_avec_canaux


def canaux_library(profile_data, width_data, height_data, thickness_data, coefficients,
                   tore_pos, debit_volumique_total_cool, nbc, plot_detail, write_in_csv, figure_dpi, plot_dir=None):
    """
    This function computes the caracteristics of channels on each point
    by interpolation between given values at injection plate (inj), end of cylindrical chamber (conv), 
    throat (col) and extremity of the nozzle (div) using library from python (scipy.interpolate)
    """

    x_value, y_value = profile_data
    lrg_inj, lrg_conv, lrg_col, lrg_tore = width_data
    ht_inj, ht_conv, ht_col, ht_tore = height_data
    e_conv, e_col, e_tore = thickness_data

    # List of x where there are channels (before the manifold) (in m)
    xcanaux = []
    # List of y where there are channels (before the manifold) (in m)
    y_coord_avec_canaux = []
    # WARNING ! ycanaux will later in this function become y on the coolant side (in m)
    i = 0
    while i < len(x_value) and x_value[i] <= tore_pos:
        xcanaux.append(x_value[i])
        y_coord_avec_canaux.append(y_value[i])
        i += 1

    pos_conv = 0  # Index of the end of cylindrical chamber
    while y_coord_avec_canaux[pos_conv] == y_coord_avec_canaux[pos_conv + 1]:
        pos_conv += 1
    pos_col = y_coord_avec_canaux.index(
        min(y_coord_avec_canaux))  # Index of the throat

    longc = len(xcanaux)  # Number of points for channels (end at the manifold)

    x_interpolate = [xcanaux[0], xcanaux[pos_conv],
                     xcanaux[pos_col], xcanaux[longc - 1]]

    y_e = [e_conv, e_conv, e_col, e_tore]
    # Thickness of the chamber wall as a function of the engine axis (in m)
    wall_thickness = [x for x in PchipInterpolator(
        x_interpolate, y_e)(xcanaux)]

    angulaire = [0]
    # y of wall on coolant side (matched with engine geometry)
    ycanaux = [y_coord_avec_canaux[0] + wall_thickness[0]]
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
        plt.plot(xcanaux, y_coord_avec_canaux,
                 color='chocolate', label='y on hot gas side')
        plt.plot(xcanaux, ycanaux, color='blue', label='y on coolant side')
        plt.title('y coordinate of wall as a function of the engine axis')
        plt.legend(loc='lower left')
        if plot_dir is not None:
            plt.savefig(f"{plot_dir}/Wall y.png")
            plt.close()
        else:
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
        if plot_dir is not None:
            plt.savefig(f"{plot_dir}/Channel thickness verif.png")
            plt.close()
        else:
            plt.show()

        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, ycanaux, color='chocolate')
        plt.title(
            'Channel travel as a function of the engine axis (y of the cold wall)')
        if plot_dir is not None:
            plt.savefig(f"{plot_dir}/Channel travel.png")
            plt.close()
        else:
            plt.show()

    debit_volumique_canal = debit_volumique_total_cool / \
                            nbc  # Volumic flow rate in a channel

    y_l = [lrg_inj, lrg_conv, lrg_col, lrg_tore]
    # Width of a channel as a function of the engine axis (in m)
    larg_canal = [x for x in PchipInterpolator(x_interpolate, y_l)(xcanaux)]
    # Width of a rib as a function of the engine axis (in m)
    larg_ailette = [(ycanaux[i] * 2 * np.pi / nbc) - larg_canal[i]
                    for i in range(0, longc)]

    y_h = [ht_inj, ht_conv, ht_col, ht_tore]
    # Height of a channel as a function of the engine axis (in m)
    ht_canal = [x for x in PchipInterpolator(x_interpolate, y_h)(xcanaux)]

    # Area of a channel as a function of the engine axis (m²)
    area_channel = []
    # Velocity of coolant in a channel as a function of the engine axis (m/s)
    vitesse_coolant = []
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
            writer.writerow(
                (1000 * xcanaux[i], 1000 * (ycanaux[i]), 1000 * (larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow(
                (1000 * xcanaux[i], 1000 * (ycanaux[i]), 1000 * (-larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow(
                (1000 * xcanaux[i], 1000 * (ycanaux[i] + ht_canal[i]), 1000 * (larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])
        writer.writerow(["StartCurve"])
        for i in range(0, longc, 3):
            writer.writerow(
                (1000 * xcanaux[i], 1000 * (ycanaux[i] + ht_canal[i]), 1000 * (- larg_canal[i] / 2)))
        writer.writerow(["EndCurve"])

        # connection between two adjacent edges (exit)
        writer.writerow(["StartCurve"])
        writer.writerow((1000 * xcanaux[0], 1000 * ycanaux[0], 500 * larg_canal[0]))
        writer.writerow((1000 * xcanaux[0], 1000 * ycanaux[0], -500 * larg_canal[0]))
        writer.writerow(["EndCurve"])

        writer.writerow(["StartCurve"])
        writer.writerow((1000 * xcanaux[0], 1000 * ycanaux[0], 500 * larg_canal[0]))
        writer.writerow((1000 * xcanaux[0], 1000 * (ycanaux[0] + ht_canal[0]), 500 * larg_canal[0]))
        writer.writerow(["EndCurve"])

        writer.writerow(["StartCurve"])
        writer.writerow((1000 * xcanaux[0], 1000 * (ycanaux[0] + ht_canal[0]), 500 * larg_canal[0]))
        writer.writerow((1000 * xcanaux[0], 1000 * (ycanaux[0] + ht_canal[0]), -500 * larg_canal[0]))
        writer.writerow(["EndCurve"])

        writer.writerow(["StartCurve"])
        writer.writerow((1000 * xcanaux[0], 1000 * (ycanaux[0] + ht_canal[0]), -500 * larg_canal[0]))
        writer.writerow((1000 * xcanaux[0], 1000 * ycanaux[0], -500 * larg_canal[0]))
        writer.writerow(["EndCurve"])

        # connection between two adjacent edges (injection)
        writer.writerow(["StartCurve"])
        writer.writerow((1000 * xcanaux[-1], 1000 * ycanaux[-1], 500 * larg_canal[-1]))
        writer.writerow((1000 * xcanaux[-1], 1000 * ycanaux[-1], -500 * larg_canal[-1]))
        writer.writerow(["EndCurve"])

        writer.writerow(["StartCurve"])
        writer.writerow((1000 * xcanaux[-1], 1000 * ycanaux[-1], 500 * larg_canal[-1]))
        writer.writerow((1000 * xcanaux[-1], 1000 * (ycanaux[-1] + ht_canal[-1]), 500 * larg_canal[-1]))
        writer.writerow(["EndCurve"])

        writer.writerow(["StartCurve"])
        writer.writerow((1000 * xcanaux[-1], 1000 * (ycanaux[-1] + ht_canal[-1]), 500 * larg_canal[-1]))
        writer.writerow((1000 * xcanaux[-1], 1000 * (ycanaux[-1] + ht_canal[-1]), -500 * larg_canal[-1]))
        writer.writerow(["EndCurve"])

        writer.writerow(["StartCurve"])
        writer.writerow((1000 * xcanaux[-1], 1000 * (ycanaux[-1] + ht_canal[-1]), -500 * larg_canal[-1]))
        writer.writerow((1000 * xcanaux[-1], 1000 * ycanaux[-1], -500 * larg_canal[-1]))
        writer.writerow(["EndCurve"])

        writer.writerow(["End"])
        file.close()

    if plot_detail >= 3:
        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, larg_ailette, label='Rib width', color='chocolate')
        plt.plot(xcanaux, larg_canal, label='Channel width', color='green')
        plt.title('Width of channels and ribs')
        plt.legend(loc='upper left')
        if plot_dir is not None:
            plt.savefig(f"{plot_dir}/Width of channels and ribs.png")
            plt.close()
        else:
            plt.show()

        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, wall_thickness, color='chocolate')
        plt.title('Thickness of chamber wall as a function of the engine axis')
        if plot_dir is not None:
            plt.savefig(f"{plot_dir}/Thickness of chamber wall.png")
            plt.close()
        else:
            plt.show()

        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, ht_canal, color='chocolate')
        plt.title('Channel height as a function of the engine axis')
        if plot_dir is not None:
            plt.savefig(f"{plot_dir}/Channel height.png")
            plt.close()
        else:
            plt.show()

    if plot_detail >= 1:
        plt.figure(dpi=figure_dpi)
        plt.plot(xcanaux, area_channel, color='chocolate')
        plt.title('Channel cross-sectionnal area as a function of the engine axis')
        if plot_dir is not None:
            plt.savefig(f"{plot_dir}/Channel area.png")
            plt.close()
        else:
            plt.show()

    return xcanaux, ycanaux, larg_canal, larg_ailette, ht_canal, wall_thickness, area_channel, longc, y_coord_avec_canaux

import casthermo as ct
import numpy as np
import matplotlib.pyplot as plt


def carto2D(larg_wall, larg_channel, ep_wall, ep_rib, dx, Hg, lamb, Tg, Hl, Tl, marker_size, display, legend_type, location, return_temp):
    npx_wall = round(larg_wall / (2 * dx) + 1)  # Number of point in the x axis at the wall
    npy_wall = round(ep_wall / dx + 1)  # Number of point in the y axis at the wall
    npx_rib = round(((larg_wall - larg_channel) / (2 * dx)) + 1)  # Number of points in the x axis at the rib
    npy_rib = round(ep_rib / dx)  # Number of point in the y axis at the rib
    nb_point = npx_wall * npy_wall + npx_rib * npy_rib  # Number of point to compute

    # Definition of configuration type, orientation, symetry and index of surrounded points
    listing = []  # List with the caracteristics useful to build the graphic for each point
    coord = []  # Coordinate of each point related to listing
    lenl = 0  # Length of the list "listing" that will be incremented in parallel with it

    # Computation of the temperature at the surface of the hot wall
    listing.append([2, 3, 1, [-1, -1, npx_wall, 1]])
    lenl += 1
    coord.append([0, 0])
    for i in range(1, npx_wall - 1):
        o = lenl - 1
        s = lenl + npx_wall
        e = lenl + 1
        listing.append([2, 3, 0, [-1, o, s, e]])
        lenl += 1
        coord.append([i * dx, 0])
    listing.append([2, 3, 2, [-1, lenl - 1, lenl + npx_wall, -1]])
    lenl += 1
    coord.append([(npx_wall - 1) * dx, 0])

    # Computation of raws in the wall
    for h in range(1, npy_wall - 1):
        listing.append([3, 0, 1, [lenl - npx_wall, -1, lenl + npx_wall, lenl + 1]])
        lenl += 1
        coord.append([0, h * dx])
        for i in range(1, npx_wall - 1):
            n = lenl - npx_wall
            o = lenl - 1
            s = lenl + npx_wall
            e = lenl + 1
            listing.append([3, 0, 0, [n, o, s, e]])
            lenl += 1
            coord.append([i * dx, h * dx])
        listing.append([3, 0, 2, [lenl - npx_wall, lenl - 1, lenl + npx_wall, -1]])
        lenl += 1
        coord.append([(npx_wall - 1) * dx, h * dx])
    
    begin_coolant = lenl - 1
    
    # Last raw of the wall
    listing.append([2, 1, 1, [lenl - npx_wall, -1, -1, lenl + 1]])
    lenl += 1
    coord.append([0, (npy_wall - 1) * dx])
    for i in range(1, npx_wall - npx_rib):
        n = lenl - npx_wall
        o = lenl - 1
        e = lenl + 1
        listing.append([2, 1, 0, [n, o, -1, e]])
        lenl += 1
        coord.append([i * dx, (npy_wall - 1) * dx])
    n = lenl - npx_wall
    o = lenl - 1
    s = lenl + npx_rib
    e = lenl + 1
    listing.append([2, 1, 0, [n, o, s, e]])
    lenl += 1
    coord.append([(npx_wall - npx_rib) * dx, (npy_wall - 1) * dx])
    for i in range(npx_wall - npx_rib + 1, npx_wall - 1):
        n = lenl - npx_wall
        o = lenl - 1
        s = lenl + npx_rib
        e = lenl + 1
        listing.append([3, 0, 0, [n, o, s, e]])
        lenl += 1
        coord.append([i * dx, (npy_wall - 1) * dx])
    listing.append([3, 0, 2, [lenl - npx_wall, lenl - 1, lenl + npx_rib, -1]])
    lenl += 1
    coord.append([(npy_wall - 1) * dx, (npy_wall - 1) * dx])

    # Computation of raws in the rib
    for h in range(1, npy_rib):
        listing.append([2, 2, 0, [lenl - npx_rib, -1, lenl + npx_rib, lenl + 1]])
        lenl += 1
        coord.append([(npx_wall - npx_rib) * dx, (h + npy_wall - 1) * dx])
        for i in range(1, npx_rib - 1):
            n = lenl - npx_rib
            o = lenl - 1
            s = lenl + npx_rib
            e = lenl + 1
            listing.append([3, 0, 0, [n, o, s, e]])
            lenl += 1
            coord.append([(i + npx_wall - npx_rib) * dx, (h + npy_wall - 1) * dx])
        listing.append([3, 0, 2, [lenl - npx_rib, lenl - 1, lenl + npx_rib, -1]])
        lenl += 1
        coord.append([(npx_rib - 1 + npx_wall - npx_rib) * dx, (h + npy_wall - 1) * dx])

    # Last raw of the rib
    listing.append([4, 0, 0, [lenl - npx_rib, -1, -1, lenl + 1]])
    lenl += 1
    coord.append([(npx_wall - npx_rib) * dx, (npy_rib + npy_wall - 1) * dx])
    for i in range(1, npx_rib - 1):
        n = lenl - npx_rib
        o = lenl - 1
        e = lenl + 1
        listing.append([5, 0, 0, [n, o, -1, e]])
        lenl += 1
        coord.append([(i + npx_wall - npx_rib) * dx, (npy_rib + npy_wall - 1) * dx])
    listing.append([5, 0, 2, [lenl - npx_rib, lenl - 1, -1, -1]])
    lenl += 1
    coord.append([(npx_rib - 1 + npx_wall - npx_rib) * dx, (npy_rib + npy_wall - 1) * dx])

    # Invertible matrix solving
    reso = np.zeros(shape=(nb_point, nb_point))
    membre = np.zeros(shape=(nb_point, 1))
    for k in range(0, lenl):
        if k < begin_coolant:
            h = Hg
            Tf = Tg
            inv = 1
        else:
            h = Hl
            Tf = Tl
            inv = 1
        # Coefficient resolution
        if listing[k][0] == 1:
            a, b, c, d, x, plus = ct.cas1(h, dx, lamb, Tf, inv)
        elif listing[k][0] == 2:
            a, b, c, d, x, plus = ct.cas2(h, dx, lamb, Tf, inv)
        elif listing[k][0] == 3:
            a, b, c, d, x, plus = ct.cas3(h, dx, lamb, Tf, inv)
        elif listing[k][0] == 4:
            a, b, c, d, x, plus = ct.cas4(h, dx, lamb, Tf, inv)
        elif listing[k][0] == 5:
            a, b, c, d, x, plus = ct.cas5(h, dx, lamb, Tf, inv)
        # Orientation resolution
        if listing[k][1] == 0:
            coef1 = a
            coef2 = b
            coef3 = c
            coef4 = d
        elif listing[k][1] == 1:
            coef1 = b
            coef2 = c
            coef3 = d
            coef4 = a
        elif listing[k][1] == 2:
            coef1 = c
            coef2 = d
            coef3 = a
            coef4 = b
        elif listing[k][1] == 3:
            coef1 = d
            coef2 = a
            coef3 = b
            coef4 = c

        # Symetry resolution
        if listing[k][2] == 1:
            coef4 += coef2
        elif listing[k][2] == 2:
            coef2 += coef4

        # Placement deduction of coefficients compared to surrounded points
        insert = []
        pos = 1
        for z in listing[k][3]:
            if 0 <= z <= (nb_point - 1):
                if pos == 1:
                    insert.append([coef1, z])
                elif pos == 2:
                    insert.append([coef2, z])
                elif pos == 3:
                    insert.append([coef3, z])
                elif pos == 4:
                    insert.append([coef4, z])
                else:
                    print("error placing")
            pos += 1

        # Introduction of coefficients in the matrix
        for values in insert:
            implacement = int(values[1])
            reso[k][implacement] = values[0]
            reso[k][k] = x
            membre[k] = plus

    reso_inv = np.linalg.inv(reso)
    T = np.dot(reso_inv, membre)

    abcisse = [m[0] for m in coord]
    ordonnee = [-m[1] for m in coord]
    temperature = [t[0] for t in T]

    if display:
        moyT_hotwall = round(sum(temperature[:npx_wall]) / npx_wall)
        moyT_coolant = round(sum(temperature[begin_coolant-1:begin_coolant-1+int(npx_wall-npx_rib)]) / (npx_wall-npx_rib))
        maxT = round(max(temperature))

        tg_avg = f"{moyT_hotwall}" if len(f'{moyT_hotwall}') == 4 else f" {moyT_hotwall}"
        tl_avg = f"{moyT_coolant}" if len(f'{moyT_coolant}') == 4 else f" {moyT_coolant}"
        t_max = f"{maxT}" if len(f'{maxT}') == 4 else f" {maxT}"
        print(f"█ Mean wall temperature at hot gaz side = {tg_avg} K                           █")
        print(f"█ Mean wall temperature at coolant side = {tl_avg} K                           █")
        print(f"█ Maximum temperature in the wall       = {t_max} K                           █")
        print("█                                                                          █")

        if legend_type == 1:
            a1 = 0.003
            a2 = -0.0025
            a3 = 0.00025
            a4 = -0.0025
            a5 = 0.002
            a6 = -0.0005
        elif legend_type == 2:
            a1 = 0.0015
            a2 = -0.002
            a3 = 0.00025
            a4 = -0.002
            a5 = 0.001
            a6 = -0.0005

        title = "2D temperature distribution (in K)" + location
        plt.figure(dpi=200)
        p = plt.scatter(abcisse, ordonnee, c=temperature, marker='s', s=marker_size, cmap='rainbow')  # rainbow#prism#flag
        plt.text(a1, a2, 'Rib', horizontalalignment='center', verticalalignment='center')
        plt.text(a3, a4, 'Coolant', horizontalalignment='center', verticalalignment='center')
        plt.text(a5, a6, 'Wall', horizontalalignment='center', verticalalignment='center')
        plt.title(title, fontsize=15)
        plt.axis("equal")
        plt.colorbar(p, shrink=0.4, aspect=15)
        plt.show()

    if return_temp:
        return [t for t in temperature[:npx_wall]]

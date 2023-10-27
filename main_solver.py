import numpy as np
import fluid_properties as flp
from tqdm import tqdm
import cte_tools as t


def mainsolver(hotgas_data, coolant_data, channel_data, chamber_data):
    """
    This is the main function used for solving the 1D case.
    The geometry is discretised into a 1 dimensionnal set of points.
    The function uses a marching algorithm, and computes all the relevant physical
    quantities at each point. The values obtained are then used on the next point.
    """

    hotgas_temp_list, molar_mass, gamma_list, Pc, c_star, PH2O, PCO2 = hotgas_data
    init_coolant_temp, init_coolant_pressure, fluid, \
    debit_mass_coolant = coolant_data
    xcanaux, ycanaux, larg_canal, larg_ailette_list, ht_canal, wall_thickness, \
    area_channel, nb_points_channel = channel_data
    y_coord_avec_canaux, nbc, diam_throat, curv_radius_pre_throat, area_throat, roughness, \
    cross_section_area_list, mach_list, material_name = chamber_data

    # Lists containing the physical quantities at each point
    coolant_temp_list = [init_coolant_temp]
    coolant_pressure_list = [init_coolant_pressure]
    coolant_viscosity_list = [flp.viscosity(init_coolant_pressure, init_coolant_temp, fluid)]
    coolant_cond_list = [flp.conductivity(init_coolant_pressure, init_coolant_temp, fluid)]
    coolant_cp_list = [flp.cp(init_coolant_pressure, init_coolant_temp, fluid)]
    coolant_density_list = [flp.density(init_coolant_pressure, init_coolant_temp, fluid)]
    hotgas_viscosity_list = []
    hotgas_cp_list = []
    hotgas_cond_list = []
    hotgas_prandtl_list = []
    coolant_reynolds_list = []
    hg_list = []
    sigma_list = []
    wall_cond_list = []
    hotwall_temp_list = []
    coldwall_temp_list = []
    flux_list = []
    coolant_velocity_list = []
    sound_speed_list = []
    hl_normal_list = []
    hl_corrected_list = []
    hl_corrected_list_2 = []
    q_list_H2O = []
    q_list_CO2 = []
    qRad_list = []
    index_throat = y_coord_avec_canaux.index(min(y_coord_avec_canaux))

    # This is to avoid oscillations near the inlet because of division by zero
    length_from_inlet = 0.03

    # Initial guess for the wall temperature
    coldwall_temp = 300
    hotwall_temp = 300

    with tqdm(total=nb_points_channel,
              desc="█ Global resolution            ",
              unit="|   █", bar_format="{l_bar}{bar}{unit}",
              ncols=76) as pbar_main:

        # Main computation loop
        for i in range(0, nb_points_channel):
            # Hydraulic diameter (4 * Area/Perimeter )
            Dhy = (2 * ht_canal[i] * larg_canal[i]) / (ht_canal[i] + larg_canal[i])

            # Velocity of the coolant
            v_cool = debit_mass_coolant / (nbc * coolant_density_list[i] * area_channel[i])
            coolant_velocity_list.append(v_cool)

            # Reynolds number of the coolant
            Re_cool = (v_cool * Dhy * coolant_density_list[i]) / coolant_viscosity_list[i]
            coolant_reynolds_list.append(Re_cool)

            # Prandtl number of the coolant
            Pr_cool = (coolant_viscosity_list[i] * coolant_cp_list[i]) / coolant_cond_list[i]

            # Compute viscosity, Cp, conductivity and Prandtl number of the hot gases
            hotgas_visc, hotgas_cp, hotgas_cond, hotgas_prandtl = t.hotgas_properties(hotgas_temp_list[i],
                                                                                      molar_mass,
                                                                                      gamma_list[i])

            # If last point in the list
            if i == nb_points_channel - 1:
                # Distance between current point and the previous (Pythagoras)
                dl = ((xcanaux[i - 1] - xcanaux[i]) ** 2 + (ycanaux[i - 1] - ycanaux[i]) ** 2) ** 0.5
            else:
                # Distance between current point and the next (Pythagoras)
                dl = ((xcanaux[i + 1] - xcanaux[i]) ** 2 + (ycanaux[i + 1] - ycanaux[i]) ** 2) ** 0.5

            length_from_inlet += dl

            # Reuse the value at previous point for a more accurate first guess (and faster convergence)
            wall_cond = 350 if i == 0 else wall_cond_list[i - 1]
            sigma = 1 if i == 0 else sigma_list[i - 1]

            # Arbitrarely create a difference to enter the loop
            new_coldwall_temp = coldwall_temp + 10
            new_hotwall_temp = hotwall_temp + 10

            # This loop's goal is to find sigma and the wall conductivity
            # It iterates until the wall temperatures have converged
            while abs(new_coldwall_temp - coldwall_temp) > 0.1 and abs(new_hotwall_temp - hotwall_temp) > 0.1:
                coldwall_temp = new_coldwall_temp
                hotwall_temp = new_hotwall_temp

                # Gas-side convective heat transfer coefficient (Bartz equation)
                hg = (0.026 / (diam_throat ** 0.2) * (((hotgas_visc ** 0.2) * hotgas_cp) / (hotgas_prandtl ** 0.6)) * (
                        (Pc / c_star) ** 0.8) * ((diam_throat / curv_radius_pre_throat) ** 0.1) * (
                              (area_throat / cross_section_area_list[i]) ** 0.9)) * sigma

                # Coolant-side convective heat transfer coefficient from Taylor (NASA TN D-4332)
                Nu = 0.023 * Re_cool ** 0.705 * Pr_cool ** 0.8 * (coldwall_temp / coolant_temp_list[i]) ** -(
                        -0.57 + 1.59 * Dhy / length_from_inlet)

                # Nusselt number correction for the channel roughness
                xi = t.darcy_weisbach(Dhy, Re_cool, roughness) / t.darcy_weisbach(Dhy, Re_cool, 0)
                roughness_correction = xi * ((1 + 1.5 * Pr_cool ** (-1 / 6) * Re_cool ** (-1 / 8) * (Pr_cool - 1)) / (
                        1 + 1.5 * Pr_cool ** (-1 / 6) * Re_cool ** (-1 / 8) * (Pr_cool * xi - 1)))

                # Compute coolant-side convective heat-transfer coefficient
                hl = Nu * roughness_correction * (coolant_cond_list[i] / Dhy)

                # Fin dimensions
                D = 2 * y_coord_avec_canaux[i]  # Diameter inside the engine
                fin_width = (np.pi * (D + ht_canal[i] + wall_thickness[i]) - nbc * larg_canal[
                    i]) / nbc  # Width of the fin

                # Correct for the fin effect (unknown source)
                m_ = ((2 * hl) / (fin_width * wall_cond)) ** 0.5
                hl_cor = hl * ((nbc * larg_canal[i]) / (np.pi * D)) + nbc * (
                        (2 * hl * wall_cond * (((np.pi * D) / nbc) - larg_canal[i])) ** 0.5) * (
                                 (np.tanh(m_ * ht_canal[i])) / (np.pi * D))

                # Correct for the fin effect (Luka Denies)
                intermediate_calc_1 = ((2 * hl * fin_width) / wall_cond) ** 0.5 * ht_canal[i] / fin_width
                nf = np.tanh(intermediate_calc_1) / intermediate_calc_1
                hl_cor2 = hl * (larg_canal[i] + 2 * nf * ht_canal[i]) / (larg_canal[i] + fin_width)

                # Compute radiative heat transfer of H2O (W) and CO2 (C) (Luka Denies)
                qW = 5.74 * ((PH2O[i] * y_coord_avec_canaux[i]) / 1e5) ** 0.3 * (hotgas_temp_list[i] / 100) ** 3.5
                qC = 4 * ((PCO2[i] * y_coord_avec_canaux[i]) / 1e5) ** 0.3 * (hotgas_temp_list[i] / 100) ** 3.5
                qRad = qW + qC

                # Computing the heat flux and wall temperatures (Luka Denies)
                flux = (hotgas_temp_list[i] - coolant_temp_list[i] + qRad / hg) / (
                        1 / hg + 1 / hl_cor + wall_thickness[i] / wall_cond)
                new_hotwall_temp = hotgas_temp_list[i] + (qRad - flux) / hg
                new_coldwall_temp = coolant_temp_list[i] + flux / hl_cor

                # Compute new value of sigma (used in the Bartz equation)
                T_hotgas_throat = hotgas_temp_list[index_throat]
                mach_hot_gases = mach_list[i]
                sigma = (((new_hotwall_temp / (2 * T_hotgas_throat)) * (
                        1 + (((gamma_list[i] - 1) / 2) * (mach_hot_gases ** 2))) + 0.5) ** -0.68) * (
                                (1 + (((gamma_list[i] - 1) / 2) * (mach_hot_gases ** 2))) ** -0.12)

                # Compute thermal conductivity of the solid at a given temperature
                wall_cond = t.conductivity(Twg=new_hotwall_temp, Twl=new_coldwall_temp, material_name=material_name)

            coldwall_temp = new_coldwall_temp
            hotwall_temp = new_hotwall_temp

            # Compute heat exchange area between two points
            # Cold-wall version (Julien)
            dA_1 = 2 * dl * (larg_canal[i] + ht_canal[i])
            # Hot-wall version (Luka Denies)
            dA_2 = (np.pi * D * dl) / nbc

            # New temperature at next point
            delta_T_coolant = ((flux * dA_1) / ((debit_mass_coolant / nbc) * coolant_cp_list[i]))
            new_coolant_temp = coolant_temp_list[i] + delta_T_coolant

            # Solving Colebrook's formula to obtain the Darcy-Weisbach friction factor
            frict_factor = t.darcy_weisbach(Dhy, Re_cool, roughness)

            # Computing pressure loss with the Darcy-Weisbach friction factor
            delta_p = 0.5 * frict_factor * (dl / Dhy) * coolant_density_list[i] * v_cool ** 2
            new_coolant_pressure = coolant_pressure_list[i] - delta_p

            # Computing the new properties of the CH4
            if new_coolant_pressure < 0:
                raise ValueError("Negative pressure ! Pressure drop is too high.")
            new_cool_visc = flp.viscosity(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)
            new_cool_cond = flp.conductivity(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)
            new_cool_cp = flp.cp(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)
            new_cool_dens = flp.density(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)
            new_cool_sound_spd = flp.sound_speed(P=new_coolant_pressure, T=new_coolant_temp, fluid=fluid)

            # Store the results
            hotgas_viscosity_list.append(hotgas_visc)
            hotgas_cp_list.append(hotgas_cp)
            hotgas_cond_list.append(hotgas_cond)
            hotgas_prandtl_list.append(hotgas_prandtl)
            hg_list.append(hg)
            hl_normal_list.append(hl)
            hl_corrected_list.append(hl_cor)
            hl_corrected_list_2.append(hl_cor2)
            hotwall_temp_list.append(hotwall_temp)
            coldwall_temp_list.append(coldwall_temp)
            flux_list.append(flux)
            sigma_list.append(sigma)
            wall_cond_list.append(wall_cond)
            coolant_pressure_list.append(new_coolant_pressure)
            coolant_temp_list.append(new_coolant_temp)
            coolant_viscosity_list.append(new_cool_visc)
            coolant_cond_list.append(new_cool_cond)
            coolant_cp_list.append(new_cool_cp)
            coolant_density_list.append(new_cool_dens)
            sound_speed_list.append(new_cool_sound_spd)
            qRad_list.append(qRad)
            q_list_CO2.append(qC)
            q_list_H2O.append(qW)
            pbar_main.update(1)

        return hl_corrected_list, hl_corrected_list_2, hotgas_viscosity_list, \
               hotgas_cp_list, hotgas_cond_list, hotgas_prandtl_list, hg_list, \
               hotwall_temp_list, coldwall_temp_list, flux_list, sigma_list, \
               coolant_reynolds_list, coolant_temp_list, coolant_viscosity_list, \
               coolant_cond_list, coolant_cp_list, coolant_density_list, \
               coolant_velocity_list, coolant_pressure_list, wall_cond_list, \
               sound_speed_list, hl_normal_list, qRad_list, q_list_CO2, q_list_H2O

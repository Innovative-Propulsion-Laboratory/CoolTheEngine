import numpy as np
import fluid_properties as flp
from tqdm import tqdm
import cte_tools as t


def solver(hotgas_data, coolant_data, channel_data, chamber_data):
    """
    This is the main function used for solving the 1D case.
    The geometry is discretised into a 1 dimensionnal set of points.
    The function uses a marching algorithm, and computes all the relevant physical
    quantities at each point. The values obtained are then used on the next point.
    """

    hotgas_recovery_temp_list, hotgas_static_temp_list, hotgas_visc_list, hotgas_pr_list, \
        hotgas_cp_list, hotgas_cond_list, MolWt, gamma_list, \
        chamber_pressure, Cstar, PH2O_list, PCO2_list = hotgas_data
    coolant_inlet_temp, coolant_inlet_pressure, coolant_name, coolant_mfr = coolant_data
    nb_channels, channel_width_list, channel_height_list, effective_fin_thickness, wall_thickness, \
        hydraulic_diameter, effective_channel_cross_section, channel_centerline, beta_list = channel_data
    z_coord_list, r_coord_list, throat_diam, throat_curv_radius, throat_area, \
        channel_roughness, cross_section_area_list, mach_list, wall_material = chamber_data

    # hotgas_recovery_temp_list = np.flip(hotgas_recovery_temp_list)
    # gamma_list = np.flip(gamma_list)
    # PH2O_list = np.flip(PH2O_list)
    # PCO2_list = np.flip(PCO2_list)
    # z_coord_list = np.flip(z_coord_list)
    # r_coord_list = np.flip(r_coord_list)
    # cross_section_area_list = np.flip(cross_section_area_list)
    # mach_list = np.flip(mach_list)
    # hotgas_visc_list = np.flip(hotgas_visc_list)
    # hotgas_pr_list = np.flip(hotgas_pr_list)
    # hotgas_cp_list = np.flip(hotgas_cp_list)
    # hotgas_cond_list = np.flip(hotgas_cond_list)

    n_points = len(z_coord_list)
    index_throat = np.argmin(np.abs(z_coord_list))

    # Compute the distance between consecutive points in the channel centerline
    dx_channels = np.diff(channel_centerline[:, 0], append=channel_centerline[:, 0][-1])
    dy_channels = np.diff(channel_centerline[:, 1], append=channel_centerline[:, 1][-1])
    dz = np.diff(channel_centerline[:, 2], append=channel_centerline[:, 2][-1])
    dr = np.diff(r_coord_list, append=r_coord_list[-1])
    dl_channel_list = np.sqrt(dx_channels**2 + dy_channels**2 + dz**2)
    dl_chamber_list = np.sqrt(dr**2 + dz**2)

    # Lists containing the physical quantities at each point
    coolant_temp_list = np.zeros(n_points)
    coolant_temp_list[-1] = coolant_inlet_temp
    coolant_pressure_list = np.zeros(n_points)
    coolant_pressure_list[-1] = coolant_inlet_pressure
    coolant_viscosity_list = np.zeros(n_points)
    coolant_viscosity_list[-1] = flp.viscosity(coolant_inlet_pressure, coolant_inlet_temp, coolant_name)
    coolant_cond_list = np.zeros(n_points)
    coolant_cond_list[-1] = flp.conductivity(coolant_inlet_pressure, coolant_inlet_temp, coolant_name)
    coolant_cp_list = np.zeros(n_points)
    coolant_cp_list[-1] = flp.cp(coolant_inlet_pressure, coolant_inlet_temp, coolant_name)
    coolant_density_list = np.zeros(n_points)
    coolant_density_list[-1] = flp.density(coolant_inlet_pressure, coolant_inlet_temp, coolant_name)
    coolant_reynolds_list = np.zeros(n_points)
    coolant_prandtl_list = np.zeros(n_points)
    hg_list = np.zeros(n_points)
    sigma_list = np.zeros(n_points)
    wall_cond_list = np.zeros(n_points)
    hotwall_temp_list = np.zeros(n_points)
    coldwall_temp_list = np.zeros(n_points)
    q_conv_list = np.zeros(n_points)
    coolant_velocity_list = np.zeros(n_points)
    hl_normal_list = np.zeros(n_points)
    hl_corrected_list = np.zeros(n_points)
    q_rad_list_H2O = np.zeros(n_points)
    q_rad_list_CO2 = np.zeros(n_points)
    q_rad_list = np.zeros(n_points)

    # This is to avoid oscillations near the inlet because of division by zero
    length_from_inlet = 0.0

    # Initial guess for the wall temperature
    coldwall_temp = 300
    hotwall_temp = 300

    with tqdm(total=n_points,
              desc="█ Global resolution            ",
              unit="|   █", bar_format="{l_bar}{bar}{unit}",
              ncols=76) as pbar_main:

        # Main computation loop (reverse direction)
        for i in range(n_points - 1, -1, -1):
            # Velocity of the coolant
            coolant_velocity_list[i] = coolant_mfr / (nb_channels * coolant_density_list[i] * effective_channel_cross_section[i])

            # Reynolds number of the coolant
            coolant_reynolds_list[i] = (coolant_velocity_list[i] * hydraulic_diameter[i] * coolant_density_list[i]) / coolant_viscosity_list[i]

            # Prandtl number of the coolant
            coolant_prandtl_list[i] = (coolant_viscosity_list[i] * coolant_cp_list[i]) / coolant_cond_list[i]

            length_from_inlet += dl_channel_list[i]

            # Reuse the value at previous point for a more accurate first guess (and faster convergence)
            wall_cond_list[i] = t.conductivity(Twg=300, Twl=300, material_name=wall_material) if i == n_points - 1 else wall_cond_list[i + 1]
            sigma = 1 if i == n_points - 1 else sigma_list[i + 1]

            # Arbitrarely create a difference to enter the loop
            new_coldwall_temp = coldwall_temp + 10
            new_hotwall_temp = hotwall_temp + 10
            # Ensure scalars for while loop
            new_coldwall_temp = float(new_coldwall_temp)
            new_hotwall_temp = float(new_hotwall_temp)

            # This loop's goal is to find sigma and the wall conductivity
            # It iterates until the wall temperatures have converged
            while abs(new_coldwall_temp - coldwall_temp) > 0.1 and abs(new_hotwall_temp - hotwall_temp) > 0.1:
                coldwall_temp = new_coldwall_temp
                hotwall_temp = new_hotwall_temp

                # Gas-side convective heat transfer coefficient (Bartz equation)
                hg = (0.026 / (throat_diam ** 0.2) * (((hotgas_visc_list[i] ** 0.2) * hotgas_cp_list[i]) / (hotgas_pr_list[i] ** 0.6)) * (
                    (chamber_pressure / Cstar) ** 0.8) * ((throat_diam / throat_curv_radius) ** 0.1) * (
                    (throat_area / cross_section_area_list[i]) ** 0.9)) * sigma

                # Compute the Darcy-Weisbach friction factor
                fD = t.darcy_weisbach(hydraulic_diameter[i], coolant_reynolds_list[i], channel_roughness)

                # Coolant-side convective heat transfer coefficient from Petukhov (https://doi.org/10.1016/S0065-2717(08)70153-9)
                k1 = 1 + 3.4*fD
                k2 = 11.7 + 1.8*coolant_prandtl_list[i]**(-1.0/3.0)
                Nu = (fD/8 * coolant_reynolds_list[i]*coolant_prandtl_list[i])/(k1*k2*(coolant_prandtl_list[i]**(2.0/3.0) - 1)*np.sqrt(fD/8))

                # Compute coolant-side convective heat-transfer coefficient
                hl = Nu * (coolant_cond_list[i] / hydraulic_diameter[i])

                # Correct for the fin effect (Luka Denies)
                intermediate_calc_1 = ((2 * hl * effective_fin_thickness[i]) / wall_cond_list[i]) ** 0.5 * channel_height_list[i] / effective_fin_thickness[i]
                nf = np.tanh(intermediate_calc_1) / intermediate_calc_1
                hl_cor = hl * (channel_width_list[i] + 2 * nf * channel_height_list[i]) / (channel_width_list[i] + effective_fin_thickness[i])

                # Compute radiative heat transfer of H2O (W) and CO2 (C) (Luka Denies)
                q_rad_H2O = 5.74 * ((PH2O_list[i] * r_coord_list[i]) / 1e5) ** 0.3 * (hotgas_static_temp_list[i] / 100) ** 3.5
                q_rad_CO2 = 4 * ((PCO2_list[i] * r_coord_list[i]) / 1e5) ** 0.3 * (hotgas_static_temp_list[i] / 100) ** 3.5
                q_rad = q_rad_H2O + q_rad_CO2

                # Computing the heat flux and wall temperatures (Luka Denies)
                q_conv = (hotgas_recovery_temp_list[i] - coolant_temp_list[i] + q_rad / hg) / (
                    1 / hg + 1 / hl_cor + wall_thickness / wall_cond_list[i])
                new_hotwall_temp = hotgas_recovery_temp_list[i] + (q_rad - q_conv) / hg
                new_coldwall_temp = coolant_temp_list[i] + q_conv / hl_cor

                # Compute new value of sigma (used in the Bartz equation)
                sigma = (((new_hotwall_temp / (2 * hotgas_recovery_temp_list[i])) *
                          (1 + (((gamma_list[i] - 1) / 2) * (mach_list[i] ** 2))
                           ) + 0.5) ** -0.68) * ((1 + (((gamma_list[i] - 1) / 2)
                                                       * (mach_list[i] ** 2))) ** -0.12)

                # Compute thermal conductivity of the solid at a given temperature
                wall_cond_list[i] = t.conductivity(Twg=new_hotwall_temp, Twl=new_coldwall_temp, material_name=wall_material)

            coldwall_temp = new_coldwall_temp
            hotwall_temp = new_hotwall_temp

            # Compute heat exchange area between two points
            dA = (np.pi * 2 * r_coord_list[i] * dl_chamber_list[i]) / (nb_channels * np.cos(np.deg2rad(beta_list[i])))

            # New temperature at previous point
            delta_T_coolant = ((q_conv * dA) / ((coolant_mfr / nb_channels) * coolant_cp_list[i]))
            new_coolant_temp = coolant_temp_list[i] + delta_T_coolant

            # Computing pressure loss with the Darcy-Weisbach friction factor
            delta_p = 0.5 * fD * (dl_channel_list[i] / hydraulic_diameter[i]) * coolant_density_list[i] * coolant_velocity_list[i] ** 2
            new_coolant_pressure = coolant_pressure_list[i] - delta_p

            # Computing the new properties of the coolant
            if new_coolant_pressure < 0:
                raise ValueError("Negative pressure ! Pressure drop is too high.")

            if i > 0:
                coolant_viscosity_list[i-1] = flp.viscosity(P=new_coolant_pressure, T=new_coolant_temp, fluid=coolant_name)
                coolant_cond_list[i-1] = flp.conductivity(P=new_coolant_pressure, T=new_coolant_temp, fluid=coolant_name)
                coolant_cp_list[i-1] = flp.cp(P=new_coolant_pressure, T=new_coolant_temp, fluid=coolant_name)
                coolant_density_list[i-1] = flp.density(P=new_coolant_pressure, T=new_coolant_temp, fluid=coolant_name)
                coolant_pressure_list[i-1] = new_coolant_pressure
                coolant_temp_list[i-1] = new_coolant_temp

            # Store the results
            hg_list[i] = hg
            hl_normal_list[i] = hl
            hl_corrected_list[i] = hl_cor
            hotwall_temp_list[i] = hotwall_temp
            coldwall_temp_list[i] = coldwall_temp
            q_conv_list[i] = q_conv
            sigma_list[i] = sigma
            q_rad_list[i] = q_rad
            q_rad_list_CO2[i] = q_rad_CO2
            q_rad_list_H2O[i] = q_rad_H2O

            pbar_main.update(1)

    return hl_corrected_list, hg_list, \
        hotwall_temp_list, coldwall_temp_list, q_conv_list, sigma_list, \
        coolant_reynolds_list, coolant_temp_list, coolant_viscosity_list, \
        coolant_cond_list, coolant_cp_list, coolant_density_list, \
        coolant_velocity_list, coolant_pressure_list, wall_cond_list, hg_list, \
        hl_normal_list, hl_corrected_list, q_rad_list, q_rad_list_CO2, q_rad_list_H2O

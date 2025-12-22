import numpy as np
from scipy.interpolate import interp1d, PchipInterpolator
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import brentq
import cea
import cte_tools as t


def mach_from_area_ratio_analytic(area_ratio, gamma, subsonic=True):
    """
    Fast Mach number calculation from area ratio and gamma using isentropic relations.
    Returns subsonic or supersonic solution.
    """
    def func(M):
        return (1/M) * ((2/(gamma+1)*(1 + (gamma-1)/2*M**2))**((gamma+1)/(2*(gamma-1)))) - area_ratio

    # Subsonic solution: M in (1e-6, 1)
    # Supersonic solution: M in (1, 10)
    if subsonic:
        return brentq(func, 1e-6, 1)
    else:
        return brentq(func, 1, 10)


def mach_list_from_area_ratios(area_ratios, gamma, throat_index):
    """
    Compute Mach number array for all area ratios, joining subsonic and supersonic branches at throat_index.
    """
    mach = np.zeros_like(area_ratios)
    for i in range(throat_index):
        mach[i] = mach_from_area_ratio_analytic(area_ratios[i], gamma, subsonic=True)
    for i in range(throat_index, len(area_ratios)):
        mach[i] = mach_from_area_ratio_analytic(area_ratios[i], gamma, subsonic=False)

    return mach


def pressure_solv(mach, gamma, stagnation_pressure):
    """
    Compute static pressure from stagnation pressure, mach number, and gamma using isentropic flow relations.
    """

    return stagnation_pressure * (1 + 0.5 * (gamma - 1) * mach ** 2) ** (-gamma / (gamma - 1))


def temperature_hotgas_solv(mach, gamma, stagnation_temperature):
    """
    Compute static temperature from stagnation temperature, mach number, and gamma using isentropic flow relations.
    """

    return stagnation_temperature * (1 + 0.5 * (gamma - 1) * mach ** 2) ** -1


def get_recovery_temperature(total_temp, gamma, mach, Pr):
    """
    Compute the recovery temperature (adiabatic wall temperature) [P. Pempie]
    """
    recovery_temp = total_temp * (
        (1 + (Pr ** 0.33) * ((gamma - 1) / 2) * mach ** 2) / (1 + ((gamma - 1) / 2) * mach ** 2))

    return recovery_temp


def wall_conductivity(Twg: float, Twl: float, material_name: str):
    """
    Compute the conductivity of the chamber wall, given temperature and material

    Twg : Hot wall temperature [K]
    Twl : Cold wall temperature [K]
    """

    T_avg = (Twg + Twl) / 2
    if material_name == "AlSi10Mg":
        return 170
    if material_name == "CuCrZr":
        # return -0.0269 * T_avg + 300
        return 300
    if material_name == "Inconel_718":
        return 0.0138 * T_avg + 5.577


def wall_cte(Twg: float, Twl: float, material_name: str):
    """
    Compute the Coefficient of Thermal Expansion, given temperature and material

    Twg : Hot wall temperature [K]
    Twl : Cold wall temperature [K]
    """

    T_avg = (Twg + Twl) / 2
    if material_name == "AlSi10Mg":
        return 27e-6
    if material_name == "CuCrZr":
        return 3e-9 * T_avg + 16.1e-6
    if material_name == "Inconel_718":
        return 7.33e-9 * T_avg + 9.07e-6


def wall_temp_limit(material_name: str):
    """
    Return the maximum acceptable hotwall temperature as a function of material
    """
    if material_name == "AlSi10Mg":
        return 473
    if material_name == "CuCrZr":
        return 800
    if material_name == "Inconel_718":
        return 1000


def wall_young_modulus(Twg: float, Twl: float, material_name: str):
    T_avg = (Twg + Twl) / 2

    if material_name == "AlSi10Mg":
        # Polynomial model valid between 20°C (293K) and 370°C (643K)
        return (-5.11e-4 * T_avg**2 + 0.346 * T_avg + 13.4) * 1e9

    if material_name == "CuCrZr":
        # Polynomial valid from 0°C (273K) to 800°C (1073K)
        return (-6e-5 * T_avg ** 2 + 1.47e-2 * T_avg + 130) * 1e9

    if material_name == "Inconel_718":
        # Polynomial valid from 20°C (293K) to 750°C (1023K)
        return (-1.59e-4 * T_avg ** 2 + 0.1015 * T_avg + 185) * 1e9


def wall_yield_strength(Twg: float, Twl: float, material_name: str):
    T_avg = (Twg + Twl) / 2

    if material_name == "AlSi10Mg":
        # Bi-linear model valid between 20°C (293K) and 370°C (643K)
        if T_avg <= 473:
            return (-0.452 * T_avg + 421) * 1e6
        elif T_avg > 473:
            return (-0.963 * T_avg + 661) * 1e6
        elif T_avg > 686:
            return 0

    if material_name == "CuCrZr":
        # Polynomial valid from 0°C (273K) to 600°C (773K)
        # Data from http://www.doi.org/10.2991/peee-15.2015.51
        # if T_avg < 890:
        #     return (-1.0e-3 * T_avg**2 + 0.58 * T_avg + 277) * 1e6
        # else:
        #     return 0

        if T_avg < 873:
            return (-1.17e-6*T_avg**2 + 5.5e-4*T_avg + 0.921) * 300e6  # Yield strength
            # return (-1.27e-6*T_avg**2 + 4.5e-4*T_avg + 0.952) * 420e6 # Ultimate strength
        else:
            return 0

    if material_name == "Inconel_718":
        # Bi-linear model valid between 18°C (295K) and 765°C (1038K)
        if T_avg <= 922:
            return (-0.25 * T_avg + 1255) * 1e6
        if T_avg > 922:
            return (-2.4 * T_avg + 3239) * 1e6


def hotgas_properties(gas_temp, molar_mass, gamma):
    """
    Computes the properties of the hot gases according to [INSERT SOURCE].
    """

    dyn_viscosity = 17.858 * (46.61 * 10 ** (-10)) * (molar_mass ** 0.5) * ((9 / 5) * gas_temp) ** 0.6
    Cp = 8314 * gamma / ((gamma - 1) * molar_mass)
    Lamb = dyn_viscosity * (Cp + (8314.5 / (4 * molar_mass)))
    Pr = 4 * gamma / (9 * gamma - 5)

    return dyn_viscosity, Cp, Lamb, Pr


def flux_equations(guess, *data):
    """
    Used by scipy.optimize.fsolve() to compute hot and cold wall temperature.
    """

    t_hot, t_cold = guess  # Initial guess
    hg, hl, t_g, t_c, wall_conductivity, wall_thickness = data

    # System of equations to solve
    f1 = hg * (t_g - t_hot) - (wall_conductivity / wall_thickness) * (t_hot - t_cold)
    f2 = hl * (t_cold - t_c) - (wall_conductivity / wall_thickness) * (t_hot - t_cold)

    return [f1, f2]


def darcy_weisbach(Dhy, Re, roughness):
    A0 = -0.79638*np.log((roughness/Dhy)/8.208 + 7.3357/Re)
    A1 = Re*(roughness/Dhy) + 9.3120665 * A0
    return ((8.128943 + A1)/(8.128943*A0 - 0.86859209*A1*np.log(A1 / (3.7099535*Re))))**2


def one_plot(x, y,
             xlabel=r'Default xlabel',
             ylabel=r'Default ylabel',
             xmin=None, xmax=None,
             ymin=None, ymax=None,
             title=r'Default title', equal_axes=False, show_grid=True,
             fmt='-k', lw=1.5, dpi=150, sci_notation=False, show=True,
             reverse=False):
    serif = {'fontname': 'DejaVu Serif'}

    margin = 0.05
    if xmin is None:
        xmin = min(x)
    if xmax is None:
        xmax = max(x)
    if ymin is None:
        ymin = min(y) - margin * min(y)
    if ymax is None:
        ymax = max(y) + margin * max(y)

    plt.rcParams["mathtext.fontset"] = 'dejavuserif'
    plt.rcParams['xtick.direction'] = 'inout'
    plt.rcParams['ytick.direction'] = 'inout'

    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    ax.minorticks_on()
    if reverse:
        x_rev = [x_ for x_ in reversed(x)]
        ax.plot(x_rev, y, fmt, linewidth=lw)
    else:
        ax.plot(x, y, fmt, linewidth=lw)
    ax.set_xlabel(xlabel, **serif)
    ax.set_ylabel(ylabel, **serif)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(title, **serif)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(9))
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    if equal_axes:
        ax.set_aspect('equal', adjustable='box')
    if show_grid:
        ax.grid(ls=':')
    if sci_notation:
        ax.ticklabel_format(style='sci', axis='y')

    if show:
        plt.show()

    return fig


def n_plots(x, y_list,
            y_label_list, colors_list,
            xlabel=r'Default xlabel',
            ylabel=r'Default ylabel',
            xmin=None, xmax=None,
            ymin=None, ymax=None,
            title=r'Default title', equal_axes=False, show_grid=True,
            fmt='-', lw=1.5, dpi=150,
            label_loc='best', sci_notation=False, show=True, reverse=False):
    serif = {'fontname': 'DejaVu Serif'}
    margin = 0.05
    if xmin is None:
        xmin = min(x)
    if xmax is None:
        xmax = max(x)
    if ymin is None:
        ymin = np.min(y_list)
        ymin = ymin - abs(margin * ymin)
    if ymax is None:
        ymax = np.max(y_list)
        ymax = ymax + abs(margin * ymax)

    plt.rcParams["mathtext.fontset"] = 'dejavuserif'
    plt.rcParams['xtick.direction'] = 'inout'
    plt.rcParams['ytick.direction'] = 'inout'

    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    ax.minorticks_on()
    for i, y in enumerate(y_list):
        if reverse:
            x_rev = [x_ for x_ in reversed(x)]
            ax.plot(x_rev, y, fmt, linewidth=lw,
                    label=y_label_list[i],
                    color=colors_list[i])
        else:
            ax.plot(x, y, fmt, linewidth=lw,
                    label=y_label_list[i],
                    color=colors_list[i])
    ax.set_xlabel(xlabel, **serif)
    ax.set_ylabel(ylabel, **serif)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(title, **serif)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(9))
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    ax.legend(loc=label_loc)
    if equal_axes:
        ax.set_aspect('equal', adjustable='box')
    if show_grid:
        ax.grid(ls=':')
    if sci_notation:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    if show:
        plt.show()

    return fig


def channel_around_engine_plot(x_center_list, y_center_list, z_coord_list,
                               r_coord_list, dpi=150, show=True):
    def set_axes_equal(ax):
        '''Set 3D plot axes to equal scale (so 1mm is the same on all axes).'''
        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()
        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)
        plot_radius = 0.5 * max([x_range, y_range, z_range])
        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    # 3D plot of channel centerline and surface of revolution
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111, projection='3d')
    # Plot channel centerline
    ax.plot(x_center_list, y_center_list, z_coord_list, color='blue', label='Channel centerline')
    # Plot surface of revolution
    theta = np.linspace(0, 2 * np.pi, 60)
    Z, Theta = np.meshgrid(np.flip(z_coord_list), theta)
    R = np.tile(np.flip(r_coord_list), (len(theta), 1))
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)
    ax.plot_surface(X, Y, Z, color='red', alpha=0.3, linewidth=0, antialiased=True)
    ax.set_xlabel('x (radial) [m]')
    ax.set_ylabel('y (radial) [m]')
    ax.set_zlabel('z (engine axis) [m]')
    ax.set_title('3D Channel Centerline and Engine Surface')
    ax.legend()
    set_axes_equal(ax)

    if show:
        plt.show()

    return fig


def compute_total_heat_output_hotwall(z_coord_list, r_coord_list, q_tot_list):
    # This assumes q_tot_list is per unit actual wall area (not projected area)
    dz = np.diff(z_coord_list, prepend=z_coord_list[0])
    dr = np.diff(r_coord_list, prepend=r_coord_list[0])
    dl = np.sqrt(dz**2 + dr**2)

    return np.sum(q_tot_list*2*np.pi*r_coord_list*dl)


def compute_1D_wall_stress(material_name: str, hotwall_thickness, z_coord_list,
                           r_coord_list, hotgas_static_pressure_list,
                           coolant_pressure_list, hotwall_temp_list,
                           coldwall_temp_list, nb_channels,
                           effective_fin_thickness_list, channel_height_list,
                           channel_width_list, total_wall_thickness):

    conductivity_list = np.zeros_like(z_coord_list)
    cte_list = np.zeros_like(z_coord_list)
    young_modulus_list = np.zeros_like(z_coord_list)
    cold_wall_thickness_list = np.zeros_like(z_coord_list)
    alpha_list = np.zeros_like(z_coord_list)
    equivalent_thickness_list = np.zeros_like(z_coord_list)
    membrane_stress_list = np.zeros_like(z_coord_list)
    bending_stress_list = np.zeros_like(z_coord_list)
    nu = 0.27
    Keff = 0.12

    for i in range(len(z_coord_list)):
        conductivity_list[i] = wall_conductivity(hotwall_temp_list[i], coldwall_temp_list[i], material_name)
        cte_list[i] = wall_cte(hotwall_temp_list[i], coldwall_temp_list[i], material_name)
        young_modulus_list[i] = wall_young_modulus(hotwall_temp_list[i], coldwall_temp_list[i], material_name)
        cold_wall_thickness_list[i] = total_wall_thickness - hotwall_thickness - channel_height_list[i]
        alpha_list[i] = effective_fin_thickness_list[i] / (effective_fin_thickness_list[i] + channel_width_list[i])
        equivalent_thickness_list[i] = hotwall_thickness + alpha_list[i] * cold_wall_thickness_list[i]
        membrane_stress_list[i] = hotgas_static_pressure_list[i] * r_coord_list[i] / equivalent_thickness_list[i]
        bending_stress_list[i] = Keff * (hotgas_static_pressure_list[i] - coolant_pressure_list[i]) * \
            (channel_width_list[i]+effective_fin_thickness_list[i])**2 / (hotwall_thickness**2)

    hoop_stress_list = membrane_stress_list + bending_stress_list  # Hoop stress at the inner side of the hot wall
    thermal_stress_list = -(young_modulus_list * cte_list * (hotwall_temp_list - coldwall_temp_list))/(2*(1-nu))  # Compressive stress is negative
    max_wall_stress_list = hoop_stress_list + thermal_stress_list

    return hoop_stress_list, thermal_stress_list, max_wall_stress_list


def compute_hotgas_flow(contour_file_path, nb_points, chamber_pressure, ox_mfr, fuel_mfr, ox_name, fuel_name):
    # Reading input data
    contour_data = np.genfromtxt(contour_file_path, delimiter=",", skip_header=1)
    z_coord_list = contour_data[0:, 0]/1000
    r_coord_list = contour_data[0:, 1]/1000

    # Reduce the number of points using scipy 1D interpolation
    # Create new x values evenly spaced between the min and max of the original x_coord_list
    x_new = np.linspace(z_coord_list[0], z_coord_list[-1], nb_points)
    # Interpolate y values at the new x positions
    interp_func = interp1d(z_coord_list, r_coord_list, kind='linear')
    y_new = interp_func(x_new)
    # Replace the original lists with the interpolated ones
    z_coord_list = x_new
    r_coord_list = y_new

    throat_radius = np.min(r_coord_list)  # Radius of the throat (in m)
    exit_radius = r_coord_list[-1]  # Radius of the exit (in m)
    throat_diam = 2 * throat_radius  # Diameter of the throat (in m)

    throat_area = np.pi * throat_radius**2  # Area of the throat (in m²)
    exit_area = np.pi * exit_radius**2  # Area of the exit (in m²)

    expansion_ratio = exit_area / throat_area  # Expansion ratio (Ae/A*)
    throat_curv_radius = 1.5 * throat_radius  # Curvature radius before the throat (in m)

    # Find the index of the throat
    i_throat = np.argmin(np.abs(r_coord_list))
    # Find the index of the beginning of the convergent section
    i_convergent = 0
    for i in range(1, len(r_coord_list)):
        if r_coord_list[i] < r_coord_list[i-1]:
            i_convergent = i
            break

    Cstar, Tc, MolWt = cea.compute_Cstar_Tc_MolWt(chamber_pressure, ox_mfr/fuel_mfr,
                                                  ox_name, fuel_name, expansion_ratio)

    hotgas_mu_chamber, hotgas_cp_chamber, hotgas_lambda_chamber, hotgas_pr_chamber, \
        hotgas_mu_throat, hotgas_cp_throat, hotgas_lambda_throat, hotgas_pr_throat, \
        hotgas_mu_exit, hotgas_cp_exit, hotgas_lambda_exit, hotgas_pr_exit = cea.get_hotgas_properties(chamber_pressure, ox_mfr/fuel_mfr,
                                                                                                       ox_name, fuel_name, expansion_ratio)

    x_chamber_throat_exit = [z_coord_list[0], z_coord_list[i_convergent],
                             z_coord_list[i_throat], z_coord_list[-1]]

    hotgas_visc_list = PchipInterpolator(x_chamber_throat_exit,
                                         [hotgas_mu_chamber, hotgas_mu_chamber,
                                          hotgas_mu_throat, hotgas_mu_exit])(z_coord_list)
    hotgas_cp_list = PchipInterpolator(x_chamber_throat_exit,
                                       [hotgas_cp_chamber, hotgas_cp_chamber,
                                        hotgas_cp_throat, hotgas_cp_exit])(z_coord_list)
    hotgas_cond_list = PchipInterpolator(x_chamber_throat_exit,
                                         [hotgas_lambda_chamber, hotgas_lambda_chamber,
                                          hotgas_lambda_throat, hotgas_lambda_exit])(z_coord_list)
    hotgas_pr_list = PchipInterpolator(x_chamber_throat_exit,
                                       [hotgas_pr_chamber, hotgas_pr_chamber,
                                        hotgas_pr_throat, hotgas_pr_exit])(z_coord_list)

    # Computation of the cross-sectional area along the engine
    cross_section_area_list = np.pi*r_coord_list**2

    gamma_list = cea.compute_gamma(chamber_pressure, ox_mfr/fuel_mfr,
                                   ox_name, fuel_name,
                                   cross_section_area_list/throat_area)

    # Computation of mach number of the hot gases
    mach_list = t.mach_list_from_area_ratios(cross_section_area_list/throat_area,
                                             gamma_list[0], i_throat)

    # Static pressure computation
    static_pressure_list = np.zeros_like(z_coord_list)  # (in Pa)
    for i in range(0, nb_points):
        static_pressure_list[i] = t.pressure_solv(mach_list[i], gamma_list[i],
                                                  chamber_pressure)

    # Partial pressure computation and interpolation of the molar fraction
    molFracH2O_chamber, molFracH2O_throat, molFracH2O_exit = cea.compute_H2O_molar_fractions(chamber_pressure, ox_mfr/fuel_mfr,
                                                                                             ox_name, fuel_name, expansion_ratio)
    molFracCO2_chamber, molFracCO2_throat, molFracCO2_exit = cea.compute_CO2_molar_fractions(chamber_pressure, ox_mfr/fuel_mfr,
                                                                                             ox_name, fuel_name, expansion_ratio)

    # Linear interpolation of molar fractions for H2O and CO2
    molFracH2O = PchipInterpolator(x_chamber_throat_exit, [molFracH2O_chamber,
                                                           molFracH2O_chamber,
                                                           molFracH2O_throat,
                                                           molFracH2O_exit])(z_coord_list)
    molFracCO2 = PchipInterpolator(x_chamber_throat_exit, [molFracCO2_chamber,
                                                           molFracCO2_chamber,
                                                           molFracCO2_throat,
                                                           molFracCO2_exit])(z_coord_list)

    # Partial pressure of H2O and CO2
    P_H2O_list = np.array([static_pressure_list[i] * molFracH2O[i] for i in range(0, nb_points)])
    P_CO2_list = np.array([static_pressure_list[i] * molFracCO2[i] for i in range(0, nb_points)])

    # Hot gas temperature computation
    hotgas_static_temp_list = np.zeros_like(z_coord_list)  # List of static hot gas temperatures
    for i in range(0, nb_points):
        hotgas_static_temp_list[i] = t.temperature_hotgas_solv(mach_list[i], gamma_list[i], Tc)

    # We assume that the total temperature is constant
    hotgas_total_temp_list = Tc*np.ones_like(z_coord_list)  # List of total hot gas temperatures

    # Computation of adiabatic wall temperature (recovery temperature)
    hotgas_recovery_temp_list = np.array([t.get_recovery_temperature(hotgas_total_temp_list[i], gamma_list[i], mach_list[i], hotgas_pr_list[i]) for i in
                                          range(0, nb_points)])

    return (z_coord_list, r_coord_list, throat_diam, throat_area,
            throat_curv_radius, Cstar, Tc, MolWt, hotgas_visc_list, hotgas_cp_list,
            hotgas_cond_list, hotgas_pr_list, cross_section_area_list, gamma_list,
            mach_list, static_pressure_list, molFracCO2, molFracH2O, P_CO2_list, P_H2O_list,
            hotgas_static_temp_list, hotgas_total_temp_list, hotgas_recovery_temp_list, x_chamber_throat_exit)

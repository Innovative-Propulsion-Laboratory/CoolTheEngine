import numpy as np
from skaero.gasdynamics import isentropic
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import brentq

# def mach_solve(area_i, area_throat, gamma, subsonic):
#     # if area_1 == area_2:
#     #     solution_mach = mach_1
#     # else:
#     #     ome = (gamma + 1) / (2 * (gamma - 1))
#     #     part_2 = (area_1 / area_2) * (mach_1 / ((1 + ((gamma - 1) / 2) * mach_1 * mach_1) ** ome))
#     #     mach_2 = mach_1
#     #     liste = []
#     #     mach = []
#     #     # Search of the mach_2 for which part_1 is minimum (750 iterations was ideal when tested)
#     #     for i in range(0, 10000):
#     #         mach_2 += 0.00001
#     #         part_1 = mach_2 * ((1 + (((gamma - 1) / 2) * mach_2 * mach_2)) ** (-ome))
#     #         liste.append(abs(part_1 - part_2))
#     #         mach.append(mach_2)
#     #     solution_mach = mach[liste.index(min(liste))]
#     # return solution_mach
#     fl = isentropic.IsentropicFlow(gamma=gamma)
#     area_ratio = area_i / area_throat
#     solution_subsonic, solution_supersonic = isentropic.mach_from_area_ratio(fl, area_ratio)
#     if subsonic:
#         return solution_subsonic
#     else:
#         return solution_supersonic


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
        return 155
    if material_name == "CuCrZr":
        return -0.0269 * T_avg + 350
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
        return 1.13e-7 * T_avg - 2.08e-5
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
        return 750
    if material_name == "Inconel_718":
        return 1000


def wall_young_modulus(Twg: float, Twl: float, material_name: str):
    T_avg = (Twg + Twl) / 2

    if material_name == "AlSi10Mg":
        # Polynomial model valid between 20°C (293K) and 370°C (643K)
        return (-5.11e-4 * T_avg**2 * 0.346 * T_avg + 13.4) * 1e9

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
        if T_avg < 890:
            return (-1.0e-3 * T_avg**2 + 0.58 * T_avg + 277) * 1e6
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
    # def colebrook(f):
    #     return 1.0 / (-2.0 * np.log10((roughness / (Dhy * 3.7)) + 2.51 / (Re * np.sqrt(f))))**2 - f

    # sol = root_scalar(colebrook, bracket=[1e-6, 0.1], method='brentq')
    # if not sol.converged:
    #     raise RuntimeError("Colebrook equation did not converge")
    # return sol.root
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
        xmin = min(x) - margin * min(x)
    if xmax is None:
        xmax = max(x) + margin * max(x)
    if ymin is None:
        ymin = np.min(y_list)
        ymin = ymin - margin * ymin
    if ymax is None:
        ymax = np.max(y_list)
        ymax = ymax + margin * ymax

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


def compute_1D_wall_stress(material_name: str, wall_thickness, z_coord_list,
                           r_coord_list, static_pressure_list,
                           coolant_pressure_list, hotwall_temp_list,
                           coldwall_temp_list):

    conductivity_list = np.zeros_like(z_coord_list)
    cte_list = np.zeros_like(z_coord_list)
    young_modulus_list = np.zeros_like(z_coord_list)
    nu = 0.33

    for i in range(len(z_coord_list)):
        conductivity_list[i] = wall_conductivity(hotwall_temp_list[i], coldwall_temp_list[i], material_name)
        cte_list[i] = wall_cte(hotwall_temp_list[i], coldwall_temp_list[i], material_name)
        young_modulus_list[i] = wall_young_modulus(hotwall_temp_list[i], coldwall_temp_list[i], material_name)

    hoop_stress_list = (static_pressure_list - coolant_pressure_list) * r_coord_list/wall_thickness
    thermal_stress_list = (young_modulus_list * cte_list * (hotwall_temp_list - coldwall_temp_list))/(2*(1-nu))
    max_wall_stress_list = np.abs(hoop_stress_list) + thermal_stress_list

    return hoop_stress_list, thermal_stress_list, max_wall_stress_list

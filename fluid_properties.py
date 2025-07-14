from CoolProp.CoolProp import PropsSI
from scipy.interpolate import RegularGridInterpolator
import numpy as np


def interpolated_Tsat(P, Tsat_1D_table):
    P_grid, Tsat_values = Tsat_1D_table
    f = RegularGridInterpolator((P_grid,), Tsat_values, bounds_error=False,
                                fill_value=None)
    return float(f([P]))


def interpolated_density(P, T, dens_2D_table):
    P_grid, T_grid, vals = dens_2D_table
    f = RegularGridInterpolator((P_grid, T_grid), vals,
                                bounds_error=False, fill_value=None)
    return float(f([P, T]))


def interpolated_cp(P, T, cp_2D_table):
    P_grid, T_grid, vals = cp_2D_table
    f = RegularGridInterpolator((P_grid, T_grid), vals,
                                bounds_error=False, fill_value=None)
    return float(f([P, T]))


def interpolated_conductivity(P, T, cond_2D_table):
    P_grid, T_grid, vals = cond_2D_table
    f = RegularGridInterpolator((P_grid, T_grid), vals,
                                bounds_error=False, fill_value=None)
    return float(f([P, T]))


def interpolated_viscosity(P, T, visc_2D_table):
    P_grid, T_grid, vals = visc_2D_table
    f = RegularGridInterpolator((P_grid, T_grid), vals,
                                bounds_error=False, fill_value=None)
    return float(f([P, T]))


def density(P, T, fluid):
    return PropsSI("D", "T", T, "P", P, fluid)


def cp(P, T, fluid):
    return PropsSI("C", "P", P, "T", T, fluid)


def conductivity(P, T, fluid):
    return PropsSI("L", "T", T, "P", P, fluid)


def viscosity(P, T, fluid):
    return PropsSI("V", "T", T, "P", P, fluid)


def sound_speed(P, T, fluid):
    return PropsSI("A", "T", T, "P", P, fluid)


def entropy(P, T, fluid):
    return PropsSI("S", "T", T, "P", P, fluid)


def pressure(S, T, fluid):
    return PropsSI('P', 'T', T, 'S', S, fluid)


def DeltaT(P, T, fluid):
    return PropsSI("T", "P", P, "Q", 0, fluid) - T


def Pcrit(fluid):
    return PropsSI("Pcrit", fluid)


def Tcrit(fluid):
    return PropsSI("Tcrit", fluid)


def molar_mass(fluid):
    return PropsSI("MOLARMASS", fluid)


def Tsat(P, fluid):
    return PropsSI("T", "P", P, "Q", 0, fluid)


def latent_heat_vap(P, fluid):
    return PropsSI('H', 'P', P, 'Q', 1, fluid) - PropsSI('H', 'P', P, 'Q', 0, fluid)


def build_fluid_property_tables(Tmin, Tmax, Pmin, Pmax, fluid):
    """
    Builds:
      * Tsat_1D_table = (P_grid, Tsat_values)
      * prop_2D_tables = {
          'density': (P_grid, T_grid, density_values),
          'cp':       (P_grid, T_grid, cp_values),
          'cond':     (P_grid, T_grid, conductivity_values),
          'visc':     (P_grid, T_grid, viscosity_values)
        }
    All pressures are in Pa, Temperatures in K.
    """

    # Grids
    T_grid = np.linspace(Tmin, Tmax, 100)
    # convert bar→Pa
    P_grid = np.linspace(Pmin * 1e5, Pmax * 1e5, 100)

    # 1D saturation‐temperature
    Tsat_values = np.array([PropsSI("T", "P", P, "Q", 0, fluid)
                            for P in P_grid])

    # initialize 2D arrays
    shape = (len(P_grid), len(T_grid))
    density_values = np.empty(shape)
    cp_values = np.empty(shape)
    conductivity_values = np.empty(shape)
    viscosity_values = np.empty(shape)

    # fill tables
    for i, P in enumerate(P_grid):
        for j, T in enumerate(T_grid):
            density_values[i, j] = PropsSI("D", "T", T, "P", P, fluid)
            cp_values[i, j] = PropsSI("C", "P", P, "T", T, fluid)
            conductivity_values[i, j] = PropsSI("L", "T", T, "P", P, fluid)
            viscosity_values[i, j] = PropsSI("V", "T", T, "P", P, fluid)

    Tsat_1D_table = (P_grid, Tsat_values)
    prop_2D_tables = {
        'density':      (P_grid, T_grid, density_values),
        'cp':           (P_grid, T_grid, cp_values),
        'conductivity': (P_grid, T_grid, conductivity_values),
        'viscosity':    (P_grid, T_grid, viscosity_values),
    }
    return Tsat_1D_table, prop_2D_tables

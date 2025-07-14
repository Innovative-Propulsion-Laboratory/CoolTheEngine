import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import RegularGridInterpolator
import time


def build_coolprop_tables(Tmin, Tmax, Pmin_bar, Pmax_bar, dT, dP_bar, fluid):
    """
    Builds:
      - Tsat_1D:   tuple (P_grid [Pa], Tsat_values)
      - prop_2D:   dict mapping prop name -> (P_grid [Pa], T_grid [K], values)
      - interp_funs: dict mapping prop name -> RegularGridInterpolator
    All arrays Ny×Nx with dimensions (len(P_grid), len(T_grid)).
    """
    # 1. Define grids in absolute units
    #    Temperatures in K
    T_grid = np.linspace(Tmin, Tmax, int(np.round((Tmax-Tmin)/dT)) + 1)
    #    Pressures in Pa (input in bar → Pa)
    P_grid = np.linspace(Pmin_bar * 1e5,
                         Pmax_bar * 1e5,
                         int(np.round((Pmax_bar-Pmin_bar)/dP_bar)) + 1)

    # 2. Pre‐allocate arrays
    shape = (len(P_grid), len(T_grid))
    dens_table = np.empty(shape)
    cp_table = np.empty(shape)
    cond_table = np.empty(shape)
    visc_table = np.empty(shape)

    # 3. Fill the 2D tables
    for i, P in enumerate(P_grid):
        for j, T in enumerate(T_grid):
            dens_table[i, j] = PropsSI("D", "T", T, "P", P, fluid)
            cp_table[i, j] = PropsSI("C", "P", P, "T", T, fluid)
            cond_table[i, j] = PropsSI("L", "T", T, "P", P, fluid)
            visc_table[i, j] = PropsSI("V", "T", T, "P", P, fluid)

    # 4. Build the 1D saturation‐temperature curve
    Tsat_vals = np.array([
        PropsSI("T", "P", P, "Q", 0, fluid)
        for P in P_grid
    ])

    # 5. Wrap interpolators (linear, allow extrapolation)
    interp_dens = RegularGridInterpolator(
        (P_grid, T_grid), dens_table,
        bounds_error=False, fill_value=None
    )
    interp_cp = RegularGridInterpolator(
        (P_grid, T_grid), cp_table,
        bounds_error=False, fill_value=None
    )
    interp_cond = RegularGridInterpolator(
        (P_grid, T_grid), cond_table,
        bounds_error=False, fill_value=None
    )
    interp_visc = RegularGridInterpolator(
        (P_grid, T_grid), visc_table,
        bounds_error=False, fill_value=None
    )
    interp_Tsat = RegularGridInterpolator(
        (P_grid,), Tsat_vals,
        bounds_error=False, fill_value=None
    )

    # 6. Package everything up
    Tsat_1D = (P_grid, Tsat_vals, interp_Tsat)
    tables_2D = {
        'density':      (P_grid, T_grid, dens_table, interp_dens),
        'cp':           (P_grid, T_grid, cp_table, interp_cp),
        'conductivity': (P_grid, T_grid, cond_table, interp_cond),
        'viscosity':    (P_grid, T_grid, visc_table, interp_visc)
    }
    return Tsat_1D, tables_2D


def interpolated_Tsat(P, fluid, Tsat_1D):
    """
    P in Pa → Tsat in K
    """
    _, _, interp = Tsat_1D
    # f expects an array of shape (n_points, ndim); here ndim=1
    return interp([P]).item()


def interpolated_density(P, T, fluid, dens_2D):
    """
    P in Pa, T in K → density in kg/m^3
    """
    _, _, _, interp = dens_2D
    return interp([P, T]).item()


def interpolated_cp(P, T, fluid, cp_2D):
    """
    P in Pa, T in K → cp in J/(kg·K)
    """
    _, _, _, interp = cp_2D
    return interp([P, T]).item()


def interpolated_conductivity(P, T, fluid, cond_2D):
    """
    P in Pa, T in K → thermal conductivity in W/(m·K)
    """
    _, _, _, interp = cond_2D
    return interp([P, T]).item()


def interpolated_viscosity(P, T, fluid, visc_2D):
    """
    P in Pa, T in K → dynamic viscosity in Pa·s
    """
    _, _, _, interp = visc_2D
    return interp([P, T]).item()


if __name__ == "__main__":
    Tsat_1D, tables_2D = build_coolprop_tables(
        Tmin=273, Tmax=500,
        Pmin_bar=1, Pmax_bar=60,
        dT=0.5, dP_bar=0.2,
        fluid="Ethanol"
    )

    reldif_rho, reldif_mu, reldif_cp, reldif_cond = [], [], [], []
    # Exact vs interpolated
    for P_test in np.linspace(1e5, 60e5):
        for T_test in np.linspace(273, 500):
            exact_rho = PropsSI("D", "T", T_test, "P", P_test, "Ethanol")
            rho_int = interpolated_density(P_test, T_test, "Ethanol", tables_2D['density'])
            reldif_rho.append(100*np.abs(exact_rho - rho_int)/exact_rho)

            exact_mu = PropsSI("V", "T", T_test, "P", P_test, "Ethanol")
            mu_int = interpolated_viscosity(P_test, T_test, "Ethanol", tables_2D['viscosity'])
            reldif_mu.append(100*np.abs(exact_mu - mu_int)/exact_mu)

            exact_cp = PropsSI("CPMASS", "T", T_test, "P", P_test, "Ethanol")
            cp_int = interpolated_cp(P_test, T_test, "Ethanol", tables_2D['cp'])
            reldif_cp.append(100*np.abs(exact_cp - cp_int)/exact_cp)

            exact_cond = PropsSI("CONDUCTIVITY", "T", T_test, "P", P_test, "Ethanol")
            cond_int = interpolated_conductivity(P_test, T_test, "Ethanol", tables_2D['conductivity'])
            reldif_cond.append(100*np.abs(exact_cond - cond_int)/exact_cond)

    start_coolprop = time.perf_counter()
    for i in range(1000):
        PropsSI("D", "T", T_test, "P", P_test, "Ethanol")
    print(f"CoolProp : {time.perf_counter() - start_coolprop}")

    start_interp = time.perf_counter()
    for i in range(1000):
        interpolated_density(P_test, T_test, "Ethanol", tables_2D['density'])
    print(f"Interpolation : {time.perf_counter() - start_interp}")

    print(np.average(reldif_rho))
    print(np.max(reldif_rho))
    print(np.average(reldif_mu))
    print(np.max(reldif_mu))
    print(np.average(reldif_cp))
    print(np.max(reldif_cp))
    print(np.average(reldif_cond))
    print(np.max(reldif_cond))

from CoolProp.CoolProp import PropsSI
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import numba as nb


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


@nb.njit
def interp1D(x, xg, yg):
    n = xg.shape[0]
    if x <= xg[0]:
        return yg[0]
    if x >= xg[-1]:
        return yg[-1]
    for i in range(n-1):
        if xg[i] <= x < xg[i+1]:
            t = (x-xg[i])/(xg[i+1]-xg[i])
            return yg[i] + t*(yg[i+1]-yg[i])
    return yg[-1]


@nb.njit
def interp2D(lp, T, lp_g, T_g, tab):
    np_, nT = lp_g.shape[0], T_g.shape[0]
    # find i
    if lp <= lp_g[0]:
        i = 0
    elif lp >= lp_g[np_-1]:
        i = np_-2
    else:
        for k in range(np_-1):
            if lp_g[k] <= lp < lp_g[k+1]:
                i = k
                break
    # find j
    if T <= T_g[0]:
        j = 0
    elif T >= T_g[nT-1]:
        j = nT-2
    else:
        for k in range(nT-1):
            if T_g[k] <= T < T_g[k+1]:
                j = k
                break

    x1, x2 = lp_g[i], lp_g[i+1]
    y1, y2 = T_g[j], T_g[j+1]
    Q11, Q21 = tab[i,  j], tab[i+1, j]
    Q12, Q22 = tab[i,  j+1], tab[i+1, j+1]
    tx = (lp-x1)/(x2-x1)
    ty = (T-y1)/(y2-y1)
    return (Q11*(1-tx)*(1-ty) + Q21*tx*(1-ty)
            + Q12*(1-tx)*ty + Q22*tx*ty)


def build_tables(Tmin, Tmax, Pmin_bar, Pmax_bar, dT, nP, fluid):
    """
    Vectorized build of:
      - T_grid [K]
      - logP_grid = log(P_grid) [for interp]
      - P_grid [Pa]
      - Tsat1D[P] → T_sat [K]
      - dens, cp, cond, visc: shape (nP, nT)
    """
    # 1) define grids
    T_grid = np.arange(Tmin, Tmax + 1e-9, dT)               # (nT,)
    P_grid = np.logspace(np.log10(Pmin_bar*1e5),
                         np.log10(Pmax_bar*1e5),
                         num=nP)                           # (nP,)
    logP_grid = np.log(P_grid)

    # 2) vectorized Tsat (will return shape (nP,))
    Tsat1D = PropsSI("T", "P", P_grid, "Q", 0, fluid)

    # 3) build the full mesh as flat arrays
    #    P_flat = [P0,P0,...,P0, P1,P1,...,P1, ..., Pn]
    #    T_flat = [T0,T1,...,Tn, T0,T1,...,Tn, ..., Tn]
    nT = T_grid.size
    P_flat = np.repeat(P_grid, nT)
    T_flat = np.tile(T_grid, nP)

    # 4) single vectorized calls for each property
    dens_flat = PropsSI("D", "P", P_flat, "T", T_flat, fluid)
    cp_flat = PropsSI("C", "P", P_flat, "T", T_flat, fluid)
    cond_flat = PropsSI("L", "T", T_flat, "P", P_flat, fluid)
    visc_flat = PropsSI("V", "T", T_flat, "P", P_flat, fluid)

    # 5) reshape back into (nP, nT) tables
    dens_table = dens_flat.reshape(nP, nT)
    cp_table = cp_flat.reshape(nP, nT)
    cond_table = cond_flat.reshape(nP, nT)
    visc_table = visc_flat.reshape(nP, nT)

    return T_grid, logP_grid, P_grid, Tsat1D, dens_table, cp_table, cond_table, visc_table


def make_interpolators(T_grid, logP_grid, P_lin, Tsat1D, dens, cp_tab, cond, visc):
    def interp_rho(P, T):
        return interp2D(np.log(P), T, logP_grid, T_grid, dens)

    def interp_mu(P, T):
        return interp2D(np.log(P), T, logP_grid, T_grid, visc)

    def interp_cp(P, T):
        return interp2D(np.log(P), T, logP_grid, T_grid, cp_tab)

    def interp_k(P, T):
        return interp2D(np.log(P), T, logP_grid, T_grid, cond)

    def interp_Tsat(P):
        return np.interp(P, P_lin, Tsat1D)

        # Prime compilation
    _ = interp_rho(5e6, 350.0)

    return interp_rho, interp_mu, interp_cp, interp_k, interp_Tsat

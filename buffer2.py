import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import numba as nb
import matplotlib.colors as mcolors
import time

# ──────────────────────────────────────────────────────────────────────────────
# 1) Build tables and Numba interpolator routines (same as before)


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


# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    fluid = "Ethanol"

    t0 = time.perf_counter()
    # 1) Build lookup tables + Numba interpolators
    T_grid, logP_grid, P_grid, Tsat1D, dens_table, cp_table, cond_table, visc_table = build_tables(
        Tmin=273, Tmax=600, Pmin_bar=1, Pmax_bar=60, dT=2, nP=100, fluid=fluid)
    print(f"Table generation time : {time.perf_counter() - t0}")

    interp_rho_func, interp_mu_func, interp_cp_func, interp_k_func, interp_Tsat = make_interpolators(
        T_grid, logP_grid, P_grid, Tsat1D, dens_table, cp_table, cond_table, visc_table)

    # Speed benchmark
    P0, T0 = 5e6, 350.0
    N = 10_000
    t0 = time.perf_counter()
    for _ in range(N):
        PropsSI("D", "T", T0, "P", P0, "Ethanol")
    t_cool = time.perf_counter() - t0

    t0 = time.perf_counter()
    for _ in range(N):
        interp_rho_func(P0, T0)
    t_nb = time.perf_counter() - t0

    print(f"CoolProp x{N}:       {t_cool:.3f} s")
    print(f"Numba-interp x{N}:    {t_nb:.3f} s   ({100*t_nb/t_cool:.1f}% of CoolProp)")

    # 2) Prepare test grid
    P_tests = np.linspace(1e5, 60e5, num=200)
    T_tests = np.linspace(273, 500, num=200)
    PP, TT = np.meshgrid(P_tests, T_tests, indexing='xy')

    # 3) Vectorize exact and approximate calls
    exact_rho_vec = np.vectorize(lambda P, T: PropsSI("D", "T", T, "P", P, fluid))
    exact_mu_vec = np.vectorize(lambda P, T: PropsSI("V", "T", T, "P", P, fluid))
    exact_cp_vec = np.vectorize(lambda P, T: PropsSI("Cpmass", "T", T, "P", P, fluid))
    exact_k_vec = np.vectorize(lambda P, T: PropsSI("conductivity", "T", T, "P", P, fluid))

    interp_rho_vec = np.vectorize(interp_rho_func)
    interp_mu_vec = np.vectorize(interp_mu_func)
    interp_cp_vec = np.vectorize(interp_cp_func)
    interp_k_vec = np.vectorize(interp_k_func)

    # 4) Compute relative-error fields in percent
    ERR = {}
    exact = exact_rho_vec(PP, TT)
    approx = interp_rho_vec(PP, TT)
    ERR['Density'] = np.abs(exact - approx) / exact * 100

    exact = exact_mu_vec(PP, TT)
    approx = interp_mu_vec(PP, TT)
    ERR['Viscosity'] = np.abs(exact - approx) / exact * 100

    exact = exact_cp_vec(PP, TT)
    approx = interp_cp_vec(PP, TT)
    ERR['Cp (mass)'] = np.abs(exact - approx) / exact * 100

    exact = exact_k_vec(PP, TT)
    approx = interp_k_vec(PP, TT)
    ERR['Conductivity'] = np.abs(exact - approx) / exact * 100

    # 5) Plot contours
    for name, err in ERR.items():
        plt.figure()
        # Mask out zeros or negatives (log scale only works for >0)
        err_pos = np.where(err > 0, err, np.nan)

        # Create a LogNorm from the smallest positive error up to the max
        vmin = np.nanmin(err_pos)
        vmax = np.nanmax(err_pos)
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)

        cs = plt.contourf(
            PP/1e5, TT, err_pos,
            levels=50,
            norm=norm,
            cmap='Spectral'
        )
        cbar = plt.colorbar(cs)
        cbar.set_label("Relative error (%) [log scale]")

        plt.title(f"{name} Relative Error (%)")
        plt.xlabel("Pressure (bar)")
        plt.ylabel("Temperature (K)")
        plt.tight_layout()

    plt.show()

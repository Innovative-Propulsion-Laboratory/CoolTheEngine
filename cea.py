from rocketcea.cea_obj import CEA_Obj
import numpy as np
from scipy.interpolate import CubicSpline


def compute_cea(Pc, MR, ox, fuel, eps):
    """
    Compute the CEA data for a given set of parameters.

    Parameters:
    Pc (float): Chamber pressure in bar.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    eps (float): Nozzle expansion ratio.

    Returns:
    dict: A dictionary containing the computed CEA data.
    """
    Pc = Pc * 14.5038  # Convert chamber pressure from bar to psi

    cea = CEA_Obj(ox=ox, fuel=fuel, Pc=Pc, MR=MR, eps=eps)
    cea.estimate_Ambient_Isp()[0]
    return {Cstar, Tc, MolWt}


def compute_gamma(Pc, MR, ox, fuel, AovAt):
    """
    Compute the adiabatic constant (gamma) for a given set of parameters.

    Parameters:
    Pc (float): Chamber pressure in psi.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    AovAt (list): Contraction/expansion ratio.

    Returns:
    list: The adiabatic constant as a function of A/At (gamma).
    """
    throat_index = np.argmin(np.abs(AovAt))  # Find the index of the throat

    cea = CEA_Obj(oxName=ox, fuelName=fuel)
    gamma = np.zeros_like(AovAt)

    # Find the first index where AovAt starts to decrease (after being constant)
    i_convergent = 0
    for i in range(1, len(AovAt)):
        if AovAt[i] < AovAt[i-1]:
            i_convergent = i
            break

    # Compute gamma[0] and gamma[throat_index] using CEA
    gamma[0] = cea.get_Chamber_MolWt_gamma(Pc=Pc, MR=MR, eps=AovAt[0])[1]
    gamma[throat_index] = cea.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=1)[1]

    # Compute gamma for points before i_convergent (if any)
    for i in range(1, i_convergent):
        gamma[i] = cea.get_Chamber_MolWt_gamma(Pc=Pc, MR=MR, eps=AovAt[i])[1]

    # Compute gamma for points after the throat
    for i, exp_ratio in enumerate(AovAt[throat_index+1:]):
        gamma[i+throat_index+1] = cea.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=exp_ratio, frozen=1, frozenAtThroat=1)[1]

    # Spline interpolation for the convergent region
    x_spline = np.arange(i_convergent, throat_index+1)
    y_spline = np.array([gamma[i_convergent-1] if i_convergent > 0 else gamma[0], gamma[throat_index]])
    cs = CubicSpline([x_spline[0], x_spline[-1]], y_spline, bc_type=((1, 0.0), 'natural'))
    gamma[i_convergent:throat_index+1] = cs(x_spline)

    return gamma


def compute_mach(Pc, MR, ox, fuel, AovAt):
    """
    Compute the Mach number for a given set of parameters.

    Parameters:
    Pc (float): Chamber pressure in psi.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    AovAt (list): Contraction/expansion ratio.

    Returns:
    list: The Mach number for every A/At.
    """
    cea = CEA_Obj(oxName=ox, fuelName=fuel)

    mach = np.zeros_like(AovAt)  # Initialize Mach number array
    throat_index = np.argmin(np.abs(AovAt))  # Find the index of the throat

    for i, cont_ratio in enumerate(AovAt[:throat_index]):
        mach[i] = cea.get_Chamber_MachNumber(Pc=Pc, MR=MR, fac_CR=cont_ratio)
    for i, exp_ratio in enumerate(AovAt[throat_index:]):
        mach[i+throat_index] = cea.get_MachNumber(Pc=Pc, MR=MR, eps=exp_ratio, frozen=1, frozenAtThroat=1)

    return mach

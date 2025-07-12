from rocketcea.cea_obj_w_units import CEA_Obj as CEA_Obj_units
from rocketcea.cea_obj import CEA_Obj
import numpy as np


def compute_Cstar_Tc_MolWt(Pc, MR, ox, fuel, exp_ratio):
    """
    Compute the CEA data for a given set of parameters.

    Parameters:
    Pc (float): Chamber pressure in bar.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    eps (float): Nozzle expansion ratio.

    Returns:
    tuple: A tuple containing the computed CEA data.
    """
    cea = CEA_Obj_units(oxName=ox, fuelName=fuel, cstar_units='m/s',
                        density_units='kg/m^3', pressure_units='Pa',
                        temperature_units='K', specific_heat_units='J/kg-K',
                        sonic_velocity_units='m/s', enthalpy_units='J/kg')

    # Compute Cstar, Tc, and MolWt using CEA methods
    Cstar = cea.get_Cstar(Pc=Pc, MR=MR)
    Tc = cea.get_Tcomb(Pc=Pc, MR=MR)
    MolWt = cea.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=exp_ratio)[0]

    return (Cstar, Tc, MolWt)


def get_hotgas_properties(Pc, MR, ox, fuel, exp_ratio):
    """
    Compute the hot gas conductivity for a given set of parameters.

    Parameters:
    Pc (float): Chamber pressure in Pa.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    exp_ratio (float): Nozzle expansion ratio.

    Returns:
    tuple: The hot gas properties
    """
    cea = CEA_Obj_units(oxName=ox, fuelName=fuel, cstar_units='m/s',
                        density_units='kg/m^3', pressure_units='Pa',
                        temperature_units='K', specific_heat_units='J/kg-K',
                        sonic_velocity_units='m/s', enthalpy_units='J/kg',
                        viscosity_units='millipoise', thermal_cond_units='W/cm-degC')

    # Get the properties of the hot gases
    cp_chamber, mu_chamber, lambda_chamber, pr_chamber = cea.get_Chamber_Transport(Pc=Pc, MR=MR, eps=exp_ratio, frozen=1)
    cp_throat, mu_throat, lambda_throat, pr_throat = cea.get_Throat_Transport(Pc=Pc, MR=MR, eps=exp_ratio, frozen=1)
    cp_exit, mu_exit,  lambda_exit, pr_exit = cea.get_Exit_Transport(Pc=Pc, MR=MR, eps=exp_ratio, frozen=1,  frozenAtThroat=1)

    # Convert units from millipoise to Pa.s for viscosity
    mu_chamber /= 10000
    mu_throat /= 10000
    mu_exit /= 10000

    # Convert thermal conductivity W/cm-K to W/m-K
    lambda_chamber *= 100
    lambda_throat *= 100
    lambda_exit *= 100

    return (mu_chamber, cp_chamber, lambda_chamber, pr_chamber,
            mu_throat, cp_throat, lambda_throat, pr_throat,
            mu_exit, cp_exit, lambda_exit, pr_exit)


def compute_gamma(Pc, MR, ox, fuel, AovAt):
    """
    Compute the adiabatic constant (gamma) for a given set of parameters.

    Parameters:
    Pc (float): Chamber pressure in Pa.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    AovAt (list): Contraction/expansion ratio.

    Returns:
    list: The adiabatic constant as a function of A/At (gamma).
    """
    cea = CEA_Obj_units(oxName=ox, fuelName=fuel, cstar_units='m/s',
                        density_units='kg/m^3', pressure_units='Pa',
                        temperature_units='K', specific_heat_units='J/kg-K',
                        sonic_velocity_units='m/s', enthalpy_units='J/kg')
    gamma = np.ones_like(AovAt)

    # Compute gamma using CEA (assumed constant along the entire engine)
    gamma = gamma*cea.get_Chamber_MolWt_gamma(Pc=Pc, MR=MR, eps=AovAt[0])[1]

    return gamma


def compute_mach(Pc, MR, ox, fuel, AovAt):
    """
    Compute the Mach number for a given set of parameters, with caching to reduce
    redundant calls to RocketCEA functions.

    Parameters:
    Pc (float): Chamber pressure in Pa.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    AovAt (list or np.ndarray): Contraction/expansion ratio.

    Returns:
    np.ndarray: The Mach number for every A/At.
    """
    Pc_psi = Pc / 1e5 * 14.5038  # Convert chamber pressure from Pa to psi
    cea = CEA_Obj(oxName=ox, fuelName=fuel)

    mach = np.zeros_like(AovAt, dtype=float)  # Initialize Mach number array
    i_t = np.argmin(np.abs(AovAt))  # Index of the throat (A/At close to 1)

    # Use dictionaries to cache previously computed results
    cr_cache = {}
    eps_cache = {}

    # Subsonic side (converging section)
    for i, cont_ratio in enumerate(AovAt[:i_t]):
        if cont_ratio not in cr_cache:
            cr_cache[cont_ratio] = cea.get_Chamber_MachNumber(Pc=Pc_psi, MR=MR, fac_CR=cont_ratio)
        mach[i] = cr_cache[cont_ratio]

    # Supersonic side (diverging section)
    for i, exp_ratio in enumerate(AovAt[i_t:]):
        if exp_ratio not in eps_cache:
            eps_cache[exp_ratio] = cea.get_MachNumber(Pc=Pc_psi, MR=MR, eps=exp_ratio, frozen=1, frozenAtThroat=1)
        mach[i + i_t] = eps_cache[exp_ratio]

    # Smooth the region near the throat
    i0 = max(i_t - 2, 0)
    i1 = min(i_t + 2, len(AovAt) - 1)

    x_old = np.array([i0, i1])
    y_old = mach[x_old]
    x_new = np.arange(i0, i1 + 1)

    # Linear interpolation for smooth transition
    y_new = np.interp(x_new, x_old, y_old)
    mach[i0:i1 + 1] = y_new

    return mach


def compute_CO2_molar_fractions(Pc, MR, ox, fuel, exp_ratio):
    """
    Compute the molar fractions of CO2 for a given set of parameters.

    Parameters:
    Pc (float): Chamber pressure in psi.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    AovAt (list): Contraction/expansion ratio.

    Returns:
    tuple: Molar fractions of CO2 at injector, throat and exit.
    """
    cea = CEA_Obj_units(oxName=ox, fuelName=fuel, cstar_units='m/s',
                        density_units='kg/m^3', pressure_units='Pa',
                        temperature_units='K', specific_heat_units='J/kg-K',
                        sonic_velocity_units='m/s', enthalpy_units='J/kg')

    molFracCO2 = cea.get_SpeciesMoleFractions(Pc=Pc, MR=MR, eps=exp_ratio)[1]['*CO2']
    return (molFracCO2[1], molFracCO2[2], molFracCO2[3])  # Return values for chamber, throat, and exit


def compute_H2O_molar_fractions(Pc, MR, ox, fuel, exp_ratio):
    """
    Compute the molar fractions of H2O for a given set of parameters.

    Parameters:
    Pc (float): Chamber pressure in psi.
    MR (float): Mixture ratio.
    ox (str): Oxidizer name.
    fuel (str): Fuel name.
    AovAt (list): Contraction/expansion ratio.

    Returns:
    tuple: Molar fractions of CO2 at injector, throat and exit.
    """
    cea = CEA_Obj_units(oxName=ox, fuelName=fuel, cstar_units='m/s',
                        density_units='kg/m^3', pressure_units='Pa',
                        temperature_units='K', specific_heat_units='J/kg-K',
                        sonic_velocity_units='m/s', enthalpy_units='J/kg')

    molFracH2O = cea.get_SpeciesMoleFractions(Pc=Pc, MR=MR, eps=exp_ratio)[1]['H2O']
    return (molFracH2O[1], molFracH2O[2], molFracH2O[3])  # Return values for chamber, throat, and exit

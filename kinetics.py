import numpy as np


def k(k0, Ea, T, R=8.3145):
    """Calculate the cinetic constant of Arrhneius.

    Args:
        k0 (float): the pre-exponencial factor.
        Ea (float): activation energy in [kJ/mol].
        T(float): Temperature in [K].
        R (float): The gas constant in [J/(mol.K)].

    Returns:
        float: The constant cinetic.
    """
    return k0 * np.exp(-Ea*1e3 / (R * T))


def fKH2S(T, k0=41769.8411, delta_Hads=2761, R=8.3145):
    """Calculate the equilibrium constant of hydrodessulfulrization.

    Args:
        T(float): Temperature in [K].
        k0 (float): the pre-exponencial factor.
        delta_Hads (float): delta energy of adsorption in [J/mol].
        R (float): The gas constant in [J/(mol.K)].

    Returns:
        float: The equilibrium constant in [cm^3/mol].

    """
    return k0 * np.exp(delta_Hads / (R * T))


def rHDS(CS, CH2, CH2S, T):
    """ Calculate the HDS reaction rate.
    Args:
        CS (float): Concentration OF sulfur compounds in [mol/cm^3].
        CH2 (float): Concentration of hydrogen in [mol/cm^3].
        CH2S (float): Concentration of hydrogen sulfide in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: The HDS reaction rate in [mol/s].
    """

    return k(4.266e9, 131.99, T) * CS * CH2**0.45/(1 + fKH2S(T) * CH2S)**2


def rHDN_B(CN_NB, CN_B, T):
    """ Calculate the HDN_B reaction rate.
    Args:
        CN_NB (float): Concentration of non-basic-nitrogen compounds in [wt%].
        CN_B (float): Concentration of basic-nitrogen in [wt%].
        T (float): Temperature in [K].

    Returns:
        float: the HDN_NB reaction rate in [wt%/s].
    """

    return  k(3.62e6, 164.94, T) * CN_NB**1.5 - k(3.66e11, 204.34, T) * CN_B**1.5 


def rHDN_NB(CN_NB, T):
    """ Calculate the HDN_NB reaction rate.
    Args:
        CN_NB (float): Concentration of non-basic-nitrogen compounds in [wt%].
        T (float): Temperature in [K].

    Returns:
        float: the HDN_B reaction rate in [wt%/s].
    """

    return k(3.62e6, 164.94, T) * CN_NB**1.5


def rHDA(pH2, CA, CNp, T):
    """ Calculate the HDA reaction rate.
    Args:
        pH2 (float): partial pressure of hydrogen in [MPa]
        CA (float): Concentration of aromatic compounds in [mol/cm^3].
        CN (float): Concentration of naphtenes compounds in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: the HDA reaction rate in [mol/s.cm^3].
    """

    return k(231.945, 80.1, T) * pH2 * CA - k(1.266e5, 112.6, T) * CNp

import numpy as np
from reactor.data import Mm1
from reactor.solubility import D1L, D2L, D4L

rhob = 0.8163 #bulk density in g/cm^3

# Component labels to use in arrays

S, NNB, NB, A, H2, H2S, Np = np.arange(7)

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


def rHDS_fun(c,T):
    """ Calculate the HDS reaction rate.
    Args:
        CS (float): Concentration of sulfur compounds in [mol/cm^3].
        CH2 (float): Concentration of hydrogen in [mol/cm^3].
        CH2S (float): Concentration of hydrogen sulfide in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: The HDS reaction rate in [mol/s.g].
    """

    return k(4.266e9, 131.99, T) * c[S] * c[H2]**0.45/(1 + fKH2S(T) * c[H2S])**2


def rHDN_B_fun(c, T, rhoL):
    """ Calculate the HDN_B reaction rate.
    Args:
        CN_NB (float): Concentration of non-basic-nitrogen compounds in [mol/cm^3].
        CN_B (float): Concentration of basic-nitrogen in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: the HDN_NB reaction rate in [mol/cm^3.s].
    """

    return  (k(3.62e6, 164.94, T) * (c[NNB] * Mm1/rhoL)**1.5 - k(3.66e11, 204.34, T) * (c[NB] * Mm1/rhoL)**1.5)
    #return (k(3.62e6, 164.94, T) * (c[NNB] * Mm1/rhoL)**1.5 - k(3.66e11, 204.34, T) * (c[NB] * Mm1/rhoL)**1.5)*rhoL/Mm1


def rHDN_NB_fun(c, T, rhoL):
    """ Calculate the HDN_NB reaction rate.
    Args:
        CN_NB (float): Concentration of non-basic-nitrogen compounds in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: the HDN_B reaction rate in [mol/cm3.s].
    """

    return (k(3.62e6, 164.94, T) * (c[NNB]* Mm1/rhoL)**1.5)
    #return  (k(3.62e6, 164.94, T) * (c[NNB]* Mm1/rhoL)**1.5)*rhoL/Mm1


def rHDA_fun(c,pH2, T, rhoL):
    """ Calculate the HDA reaction rate.
    Args:
        pH2 (float): partial pressure of hydrogen in [MPa]
        CA (float): Concentration of aromatic compounds in [mol/cm^3].
        CN (float): Concentration of naphtenes compounds in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: the HDA reaction rate in [mol/s.cm^3].
    """

    return k(231.945, 80.1, T) * pH2 * c[A] - k(1.266e5, 112.6, T) * c[Np]
    #return k(1.041e5, 121.40, T) * pH2 * c[A]*Mm1/rhoL - k(8.805e9, 186.40, T) * (1 - Mm1/rhoL*c[A])
    """Kmono = 0.7

    kdir_line = 6.04e2*np.exp(-12140/T)*10*0.85

    k_dir = kdir_line/4

    k_inv = k_dir/Kmono

    return k_dir*pH2*c[A] - k_inv*c[Np]"""






def ft_reactants(r, c, pH2, T, rhoL):

    diff_oil = D1L(T)
    diff_H2 = D2L(T)
    diff_H2S = D4L(T)

    rHDS = rHDS_fun(c,T)
    rHDN_NB = rHDN_NB_fun(c, T, rhoL)
    rHDN_B = rHDN_B_fun(c, T, rhoL)
    rHDA = rHDA_fun(c, pH2, T, rhoL)


    fS = 1 * rHDS * rhob / diff_oil
    fNB =  1 * rHDN_NB / diff_oil
    fNNB = -1 * rHDN_B / diff_oil
    fA = 1 * rHDA / diff_oil
    fH2 = (15 * rHDS * rhob + 6 * rHDN_NB + 2 * rHDN_B + 3 * rHDA) / diff_H2
    fH2S = -9 * rHDS * rhob / diff_H2S
    fNp = -1 * rHDA/ diff_oil

    return np.array([fS, fNB, fNNB, fA, fH2, fH2S, fNp])


def effective_reactions(r, c, pH2, T, rhoL):

    rr1 = rHDS_fun(c,T)
    rr2 = rHDN_NB_fun(c, T, rhoL)
    rr3 = rHDN_B_fun(c, T, rhoL)
    rr4 = rHDA_fun(c, pH2, T, rhoL)

    return np.array([rr1, rr2, rr3, rr4])

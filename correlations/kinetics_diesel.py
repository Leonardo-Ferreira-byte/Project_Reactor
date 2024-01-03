import numpy as np
from correlations.mass_transfer import get_diffusion_oil, get_diffusion_H2, get_diffusion_H2S, get_rg, get_Dk, get_De

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


def rHDS_fun(c, T):
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


def rHDN_B_fun(c, T, rhoL, Mm):
    """ Calculate the HDN_B reaction rate.
    Args:
        CN_NB (float): Concentration of non-basic-nitrogen compounds in [mol/cm^3].
        CN_B (float): Concentration of basic-nitrogen in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: the HDN_NB reaction rate in [mol/cm^3.s].
    """

    return (k(3.62e6, 164.94, T) * (c[NNB] * Mm/rhoL)**1.5 - k(3.66e11, 204.34, T) * (c[NB] * Mm/rhoL)**1.5)
    # return (k(3.62e6, 164.94, T) * (c[NNB] * Mm/rhoL)**1.5 - k(3.66e11, 204.34, T) * (c[NB] * Mm/rhoL)**1.5)*rhoL/Mm


def rHDN_NB_fun(c, T, rhoL, Mm):
    """ Calculate the HDN_NB reaction rate.
    Args:
        CN_NB (float): Concentration of non-basic-nitrogen compounds in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: the HDN_B reaction rate in [mol/cm3.s].
    """

    return (k(3.62e6, 164.94, T) * (c[NNB] * Mm/rhoL)**1.5)
    # return  (k(3.62e6, 164.94, T) * (c[NNB]* Mm/rhoL)**1.5)*rhoL/Mm


def rHDA_fun(c, pH2, T, rhoL):
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
    # return k(1.041e5, 121.40, T) * pH2 * c[A]*Mm/rhoL - k(8.805e9, 186.40, T) * (1 - Mm/rhoL*c[A])
    #return k(1.041e5, 121.40, T) * pH2 * c[A] - k(8.805e9, 186.40, T) * (c[Np])



def ft_reactants(r, c, pH2, T, rhoL, rhob, viscosity, vL, vH2, vH2S, Vg, Sg, Mm, teta):

    diff_oil = get_diffusion_oil(T, vL, viscosity)
    diff_H2 = get_diffusion_H2(T, viscosity, vL, vH2)
    diff_H2S = get_diffusion_H2S(T, viscosity, vL, vH2S)

    rg = get_rg(Vg, Sg)
    Dk = get_Dk(rg, T, Mm)
    De_oil = get_De(diff_oil, Dk, teta)
    De_H2 = get_De(diff_H2, Dk, teta)
    De_H2S = get_De(diff_H2S, Dk, teta)



    rHDS = rHDS_fun(c, T)
    rHDN_NB = rHDN_NB_fun(c, T, rhoL, Mm)
    rHDN_B = rHDN_B_fun(c, T, rhoL, Mm)
    rHDA = rHDA_fun(c, pH2, T, rhoL)

    fS = 1 * rHDS * rhob / De_oil
    fNB = 1 * rHDN_NB / De_oil
    fNNB = -1 * rHDN_B / De_oil
    fA = 1 * rHDA / De_oil
    fH2 = (15 * rHDS * rhob + 6 * rHDN_NB + 6 * rHDN_B + 3 * rHDA) / De_H2
    fH2S = -9 * rHDS * rhob / De_H2S
    fNp = -1 * rHDA / De_oil

    return np.array([fS, fNB, fNNB, fA, fH2, fH2S, fNp])


def effective_reactions(r, c, pH2, T, rhoL, Mm):

    rr1 = rHDS_fun(c, T)
    rr2 = rHDN_NB_fun(c, T, rhoL, Mm)
    rr3 = rHDN_B_fun(c, T, rhoL, Mm)
    rr4 = rHDA_fun(c, pH2, T, rhoL)

    return np.array([rr1, rr2, rr3, rr4])

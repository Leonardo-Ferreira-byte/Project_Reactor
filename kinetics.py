import numpy as np

#[0] = S

#[1] = NNB

#[2] = NB

#[3] = A

#[4] = H2

#[5] = H2S

#[6] = Np

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


def rHDS(c,T):
    """ Calculate the HDS reaction rate.
    Args:
        CS (float): Concentration of sulfur compounds in [mol/cm^3].
        CH2 (float): Concentration of hydrogen in [mol/cm^3].
        CH2S (float): Concentration of hydrogen sulfide in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: The HDS reaction rate in [mol/s].
    """

    return k(4.266e9, 131.99, T) * c[S] * c[H2]**0.45/(1 + fKH2S(T) * c[H2S])**2


def rHDN_B(c, T):
    """ Calculate the HDN_B reaction rate.
    Args:
        CN_NB (float): Concentration of non-basic-nitrogen compounds in [wt%].
        CN_B (float): Concentration of basic-nitrogen in [wt%].
        T (float): Temperature in [K].

    Returns:
        float: the HDN_NB reaction rate in [wt%/s].
    """

    #return  k(3.62e6, 164.94, T) * c[NNB]**1.5 - k(3.66e11, 204.34, T) * c[NB]**1.5 
    return 1.20/3600 * c[NNB]**1.5 -  4.85/3600 * c[NB]**1.5


def rHDN_NB(c, T):
    """ Calculate the HDN_NB reaction rate.
    Args:
        CN_NB (float): Concentration of non-basic-nitrogen compounds in [wt%].
        T (float): Temperature in [K].

    Returns:
        float: the HDN_B reaction rate in [wt%/s].
    """

    #return k(3.62e6, 164.94, T) * c[NNB]**1.5
    return  1.20/3600* c[NNB]**1.5


def rHDA(c,pH2, T):
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


def ft_reactants(r, c, pH2, T):

    ni = [-1, -1, 1, -1, -15, 9, 1]

    #checar sinal

    fS = -ni[0] * rHDS(c,T)
    fNB = -ni[1] * rHDN_NB(c, T)
    fNNB = -ni[2] * rHDN_B(c, T)
    fA = -ni[3] * rHDA(c, pH2, T)
    fH2 = -ni[4] * rHDS(c, T)
    fH2S = -ni[5] * rHDS(c, T)
    fNp = -ni[6] * rHDA(c, pH2, T)

    return np.array([fS, fNB, fNNB, fA, fH2, fH2S, fNp])


def effective_reactions(r, c, pH2, T):

    rr1 = rHDS(c,T)
    rr2 = rHDN_NB(c, T)
    rr3 = rHDN_B(c, T)
    rr4 = rHDA(c, pH2, T)

    return np.array([rr1, rr2, rr3, rr4])

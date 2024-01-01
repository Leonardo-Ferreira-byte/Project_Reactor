import numpy as np
    
def fKH2S(T, k0 = 41769.8411, delta_Hads=2761 ,R = 8.3145):
    """Calculate the equilibrium constant of hydrodessulfulrization.

    Args:
        T(float): Temperature in [K].
        k0 (float): the pre-exponencial factor.
        delta_Hads (float): delta energy of adsorption in [J/mol].
        R (float): The gas constant in [J/(mol.K)].

    Returns:
        float: The equilibrium constant in [cm^3/mol].

    """
    return  k0* np.exp(delta_Hads / (R * T))

def k_app():

    A = 0.21
    B = 1.40
    G = 0.0572
    kin = 0.67

    return (A/G**B + 1/kin)**-1

def rHDS(CS, CH2, CH2S, KH2S):
    """ Calculate the HDS reaction rate.
    Args:
        CS (float): Concentration of sulfur compounds in [mol/cm^3].
        CH2 (float): Concentration of hydrogen in [mol/cm^3].
        CH2S (float): Concentration of hydrogen sulfide in [mol/cm^3].
        T (float): Temperature in [K].

    Returns:
        float: The HDS reaction rate in [mol/s].
    """

    return k_app() * CS * CH2**0.45/(1 + KH2S * CH2S)**2
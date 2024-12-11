Patm = 0.101325  # Pressure in CNTP conditions in MPa

T0 = 273.15  # Temperature in CNTP conditions in K

alfa1 = 7  # first coeficient of gas-liquid mass transfer in cm^(-1.6).

alfa2 = 0.4  # second coeficient of gas-liquid mass transfer

roW = 0.9976  # density of water at 15.6°C in g/cm^3.

vc_2 = 65.1  # critical specific volume of hydrogen in cm^3/mol.

vc_4 = 98.6  # critical specific volume of hydrogen sulfite in cm^3/mol.

vH2 = 0.285 * vc_2**1.048  # volume molar of hydrogen in cm^3/mol.

vH2S = 0.285 * vc_4**1.048  # volume molar of hydrogen sulfite in cm^3/mol.

Vn = 8.3145 * 273.15 / 101325 * 1000  # Volume molar at standard conditions in Nl.


def specific_gravity(API):
    """Get the specific gravity of oil.
    Parameters
    ----------
    API: int or float

    Returns
    -------
    The specific_gravity of oil.
    specific_gravity: float

    """
    return 141.5 / (API + 131.5)


def get_API(specific_gravity):
    """Get the API gravity of oil.

     Parameters
     ----------
     specific_gravity: float

    Returns
    -------
    The API gravity of oil.
    API: float

    """
    return 141.5 / specific_gravity - 131.5


def get_volume_molar_oil(T_MeABP, d15_6, Mm):

    vc1_m = 7.5214e-3 * T_MeABP**0.2896 * d15_6 ** (-0.7666)
    vc1 = vc1_m * Mm * 62.4279

    return 0.285 * vc1**1.048

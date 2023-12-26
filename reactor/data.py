API = 22

alfa1 = 7 #first coeficient of gas-liquid mass transfer in cm^(-1.6).

alfa2 = 0.4 #second coeficient of gas-liquid mass transfer

T_MeABP = (476 + 273.15)*1.8 # Mean average boiling point of oil in °R.

roW = 0.9976 #density of water at 15.6°C in g/cm^3.

def specific_gravity(API):

    '''Get the specific gravity of oil.
    Parameters
    ----------
    API: int or float

    Returns
    -------
    The specific_gravity of oil.
    specific_gravity: float

    '''
    return 141.5/(API + 131.5)


d15_6 = specific_gravity(API)

Mm1 = 441.9 #molecuar weigth of the oil.

vc1_m = 7.5214e-3*T_MeABP**0.2896*d15_6**(-0.7666) #critical specific volume in ft^3/lb.

vc_2 = 65.1 #critical specific volume of hydrogen in cm^3/mol.

vc_4 = 98.6 #critical specific volume of hydrogen sulfite in cm^3/mol.

vc1 = vc1_m*Mm1*62.4279 #critical specific volume of the oil in cm^3/mol.

v1 = 0.285*vc1**1.048 #volume molar of the oil in cm^3/mol.

v2 = 0.285*vc_2**1.048   #volume molar of hydrogen in cm^3/mol.

v4 = 0.285*vc_4**1.048   #volume molar of hydrogen sulfite in cm^3/mol.

Vn = 8.3145*273.15/101325*1000 #Volume molar at standard conditions in Nl.
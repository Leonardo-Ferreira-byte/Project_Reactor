from reactor.data import Patm
def deltarhoP(rho0, P):
    '''Function to get the density of oil corretion with pressure.

        Parameters
        ----------
        rho0: int or float
             Density at standard conditions (15,6°C; 101.3 kPa) in lb/ft^3.

        P: int or float
           Pressure in psia.

        Returns
        -------
        deltarhoP: float
                  Density correction with pressure in lb/ft^3.

    '''

    return (0.167 + (16.181*10**(-0.0425*rho0))) * (P/1000) - 0.01*(0.299 + (263*10**(-0.0603*rho0))) * (P/1000)**2


def deltarhoT(rho0, P, T):
    '''Function to get the density of oil corretion with temperature.

        Parameters
        ----------
        rho0: int or float
             Density at standard conditions (15,6°C; 101.3 kPa) in lb/ft^3.

        P: int or float
           Pressure in psia.

        T: int or float
           Temperature in  °R.

        Returns
        -------
        deltarhoT: float
                  Density correction with temperature in lb/ft^3.

    '''

    return (0.0133 + 152.4*(rho0 + deltarhoP(rho0, P))**(-2.45)) * (T-520)\
    - (8.1e-6 - 0.0622*10**(-0.764*(rho0 + deltarhoP(rho0, P)))) * (T-520)**2


def oil_density(rho0, P, T):
    '''Get the density of oil in function of pressure and temperature.

        Parameters
        ----------
        rho0: int or float
             Density at standard conditions (15,6°C; 101.3 kPa) in g/cm^3.

        P: int or float
           Pressure in Pa.

        T: int or float
           Temperature in K.

        Returns
        -------
        oil_density: float
                  Density correction with temperature and pressure in g/cm^3.

    '''

    rho0 = rho0/0.016018
    T = 1.8*T
    P = P*14.695/Patm

    return (rho0 + deltarhoP(rho0, P) - deltarhoT(rho0, P, T))*0.016018

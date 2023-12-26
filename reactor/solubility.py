import numpy as np
from reactor.density_correlations import oil_density
from reactor.data import Vn, v1, v2, v4

def mi_L(T, API):

      '''Function to get the viscosity of oil.

      Parameters
      ----------
      API: float or int
            API gravity of oil.   

      T: int or float
            Temperature in K.

      Returns
      -------
      mi_L: float
            Viscosity of oil sulfite in g/(cm*s).'''

      T = T*1.8

      a = 10.313*np.log10(T - 460) - 36.447
      return 3.141e10*(T - 460)**(-3.444)*((np.log10(API))**a)/100


def D1L(T, API=22):

      '''Gets the coefficient of diffusivity of the organic sulfur compound in solution in cm^2/s.

      Parameters
      ----------
      T: int or float
            Temperature in K.

      Returns
      -------
      D1L: float
            The difusivity coefficient.
      '''

      return (8.93e-8*(v1**0.267)*T)/((v1**0.433)*(mi_L(T, API)*100))

def D2L(T, API=22):

      '''Gets the coefficient of diffusivity of hydrogen in solution in cm^2/s.

      Parameters
      ----------
      T: int or float
            Temperature in K.

      Returns
      -------
      D2L: float
            The difusivity coefficient.
      '''

      return (8.93e-8*(v1**0.267)*T)/((v2**0.433)*(mi_L(T, API)*100))

def D4L(T, API=22):

      '''Gets the coefficient of diffusivity of hydrogen sulfite in solution in cm^2/s.

      Parameters
      ----------
      T: int or float
            Temperature in K.

      Returns
      -------
      D2L: float
            The difusivity coefficient.
    '''

      return 8.93e-8*v1**0.267*T/v4**0.433/(mi_L(T, API)*100)

def Lambda2(rho0, T):

    '''Function to get the solubility of hydrogen in hydrocarbon mixtures.

    Parameters
    ----------
    T: int or float
       Temperature in K.

    Returns
    -------
    Lambda2: float
             Solubility of hydrogen in (Nl H2,)/[(g oil)*(MPa)].
    '''
    T = T - 273.15

    ro_20 = oil_density(rho0,0.101325,293.15)

    return  (-0.559729 - 0.42947e-3*T + 3.07539e-3*T/ro_20 + 1.94593e-6*T**2 + 0.835783/ro_20**2) 

def Henry_coefficient2_fun(rho0, P, T):

    '''Function to get the Henry coefficient of hydrogen in hydrocarbon mixtures.
     
        rho0: int or float
            Density at standard conditions (15,6°C; 101.3 kPa) in g/cm^3.
        
        P: int or float
           Pressure in Pa.
    
        T: int or float
           Temperature in K.

        Returns
        -------
        H: float
           The Henry coefficient for hydrogen in MPa.cm^3/mol. '''  

    
    return Vn/(Lambda2(rho0,T)*oil_density(rho0,P,T)/1000)

def Lambda4(T):

   '''Function to get the solubility of hydrogen sulfite in hidrocarbon mixtures.

    Parameters
    ----------
    T: int or float
       Temperature in °C.

    Returns
    -------
    Lambda2: float
             Solubility of hydrogen sulfite in (Nl H2S,)/[(g oil)*(MPa)].
    '''
   T = T - 273.15

   return np.exp(3.3670 - 0.008470*T)

def Henry_coefficient4_fun(rho0, P, T):

    '''Function to get the Henry coefficient of hydrogen sulfite in hydrocarbon mixtures. 

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
        H: float
           The Henry coefficient for hydrogen sulfite in MPa.cm^3/mol. '''  

    return Vn/(Lambda4(T)*oil_density(rho0,P,T)/1000)
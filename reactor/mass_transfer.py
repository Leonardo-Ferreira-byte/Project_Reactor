from reactor.data import alfa1, alfa2
from reactor.density_correlations import oil_density
from reactor.solubility import mi_L, D1L, D2L, D4L

def k2L_aL_fun(rho0, P, T, API, GL): 
    
   '''Get the gas-liquid mass transfer coefficient from liquid to particle for hydrogen.
    
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
        k2Lat: float
             Mass transfer coefficient in s^-1.'''
   
     
   return D2L(T, API)*alfa1*(GL/mi_L(T, API))**alfa2*(mi_L(T, API)/(oil_density(rho0,P,T)*D2L(T, API)))**(0.5)

def k4L_aL_fun(rho0, P, T, API, GL): 
    
   '''Get the gas-liquid mass transfer coefficient from liquid to particle for hydrogen sulfite.
    
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
        k4Lat: float
             Mass transfer coefficient in s^-1.'''
   
     
   return D4L(T, API)*alfa1*(GL/mi_L(T, API))**alfa2*(mi_L(T, API)/(oil_density(rho0,P,T)*D4L(T, API)))**(0.5)

def k1SaS_fun(rho0, P, T, API, GL, aS):

    '''Get the mass transfer coefficient from liquid to particle for organic sulfur compound.
    
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
        k1S: float
             Mass transfer coefficient from liquid to particle in s^-1.'''

    return D1L(T, API)*aS*1.8*(GL/(mi_L(T, API)*aS))**0.5*(mi_L(T, API)/(oil_density(rho0,P,T)*D1L(T, API)))**(1/3) * aS

def k2SaS_fun(rho0, P, T, API, GL, aS):

    '''Get the mass transfer coefficient from liquid to particle for hydrogen.
    
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
        k2S: float
             Mass transfer coefficient from liquid to particle in s^-1.'''

    return D2L(T, API)*aS*1.8*(GL/(mi_L(T, API)*aS))**0.5*(mi_L(T, API)/(oil_density(rho0,P,T)*D2L(T, API)))**(1/3) * aS

def k4SaS_fun(rho0, P, T, API, GL, aS):

    '''Get the mass transfer coefficient from liquid to particle for hydrogen sulfite.
    
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
        k4sat: float
             Mass transfer coefficient from liquid to particle in s^-1.'''

    return D4L(T, API) * aS*1.8 * (GL/(mi_L(T, API)*aS))**0.5 * (mi_L(T, API)/(oil_density(rho0,P,T)*D4L(T, API)))**(1/3) * aS


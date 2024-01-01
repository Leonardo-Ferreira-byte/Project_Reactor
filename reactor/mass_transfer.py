from reactor.data import alfa1, alfa2
from reactor.density_correlations import oil_density

def aS_fun(porosity, dp):

   return 6*(1-porosity)/dp

def kiL_aL_fun(Di, viscosity, rho0, P, T, GL): 
       
   return Di * alfa1 * (GL / viscosity)**alfa2 * (viscosity / (oil_density(rho0,P,T) * Di))**(0.5)

def kiS_aS_fun(Di, viscosity, rho0, P, T, GL, aS):

    return Di * aS * 1.8 * (GL / (viscosity * aS))**0.5 * (viscosity / (oil_density(rho0,P,T) * Di))**(1/3) * aS


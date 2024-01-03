from correlations.data import alfa1, alfa2
from correlations.density_correlations import oil_density

def get_diffusion_oil(T, viscosity, vL):

      return 8.93e-8 * vL**0.267 *T / vL**0.433 / viscosity / 100

def get_diffusion_H2(T, viscosity, vL, vH2):

      return 8.93e-8 * vL**0.267 * T / vH2**0.433 / viscosity / 100

def get_diffusion_H2S(T, viscosity, vL, vH2S):

      return 8.93e-8 * vL**0.267 * T/ vH2S**0.433 / viscosity / 100

def get_rg(Vg, Sg):
     
     return 2 * Vg / Sg

def get_Dk(rg, T, Mm):
     
     return 9700 * rg * (T / Mm)**0.5

def get_De(Doil, Dk, teta, tau = 0.4):

   return (teta / tau) * (1 / (1 / Doil) + 1 / (1 / Dk))

def aS_fun(porosity, dp):

   return 6*(1-porosity)/dp

def kiL_aL_fun(Di, viscosity, rho0, P, T, GL): 
       
   return Di * alfa1 * (GL / viscosity)**alfa2 * (viscosity / (oil_density(rho0,P,T) * Di))**(0.5)

def kiS_aS_fun(Di, viscosity, rho0, P, T, GL, aS):

    return Di * aS * 1.8 * (GL / (viscosity * aS))**0.5 * (viscosity / (oil_density(rho0,P,T) * Di))**(1/3) * aS


from system.reactor import Reactor
import numpy as np
from correlations.density_correlations import oil_density
from correlations.data import Patm, T0, roW, vH2, vH2S, specific_gravity, get_volume_molar_oil
from correlations.solubility import get_viscosity, get_Henry_H2, get_Henry_H2S
from correlations.mass_transfer import  get_diffusion_oil, get_diffusion_H2, get_diffusion_H2S, kiL_aL_fun, kiS_aS_fun


class Fluid:

    def __init__(self, T=653.15, P=5.3, fi=356, LHSV=2, API=22, Mm=441.9, T_MeABP=(476 + 273.15)*1.8, rhob=0.8163, heterogeneous=False):

        self.reactor = Reactor()
        self.T = T
        self.P = P
        self.fi = fi
        self.LHSV = LHSV
        self.d15_6 = specific_gravity(API)
        self.rho0 = self.d15_6 * roW
        self.rhoL = oil_density(self.rho0, P, T)
        self.Mm = Mm
        self.rhob = rhob
        self.rhop = rhob / (1 - self.reactor.porosity)
        self.teta = self.reactor.Vg * self.rhop
        self.heterogeneous = heterogeneous
        self.vH2 = vH2
        self.vH2S = vH2S
        self.vL = get_volume_molar_oil(T_MeABP, self.d15_6, self.Mm)
        self.viscosity = get_viscosity(T, API)
        self.diffusion_oil = get_diffusion_oil(T, self.viscosity, self.vL)
        self.diffusion_H2 = get_diffusion_H2(T, self.viscosity, self.vL, self.vH2)
        self.diffusion_H2S = get_diffusion_H2S(T, self.viscosity, self.vL, self.vH2S)
        self.Henry_H2 = get_Henry_H2(self.rho0, P, T)
        self.Henry_H2S = get_Henry_H2S(self.rho0, P, T)
        self.uL = self.LHSV * self.reactor.z / 3600
        self.GL = self.rhoL * self.uL
        self.uG = self.uL * (Patm / self.P) * (self.T/T0) * self.fi
        self.kLaL_H2 = kiL_aL_fun(self.diffusion_H2, self.viscosity, self.rho0, P, T, self.GL)
        self.kLaL_H2S = kiL_aL_fun(self.diffusion_H2S, self.viscosity, self.rho0, P, T, self.GL)
        self.kSaS_oil = kiS_aS_fun(self.diffusion_oil, self.viscosity, self.rho0, P, T, self.GL, self.reactor.aS)
        self.kSaS_H2 = kiS_aS_fun(self.diffusion_H2, self.viscosity, self.rho0, P, T, self.GL, self.reactor.aS)
        self.kSaS_H2S = kiS_aS_fun(self.diffusion_H2S, self.viscosity, self.rho0, P, T, self.GL, self.reactor.aS)

    def mass_balance_gas_phase_H2(self, p2G, C2L, R=8.3145):

        return - self.kLaL_H2 * (p2G / self.Henry_H2 - C2L) * R * self.T / self.uG

    def mass_balance_gas_phase_H2S(self, p4G, C4L, R=8.3145):

        return - self.kLaL_H2S * (p4G / self.Henry_H2S - C4L) * R * self.T / self.uG

    def mass_balance_gas_liquid_phase_H2(self, p2G, C2L, C2S):

        return (self.kLaL_H2 * (p2G / self.Henry_H2 - C2L) - self.kSaS_H2 * (C2L - C2S)) / self.uL

    def mass_balance_gas_liquid_phase_H2S(self, p4G, C4L, C4S):

        return (self.kLaL_H2S * (p4G / self.Henry_H2S - C4L) - self.kSaS_H2S * (C4L - C4S)) / self.uL

    def mass_balance_liquid(self, C1L, C1S):

        return - (C1L - C1S) * self.kSaS_oil / self.uL

    def wt_to_molar(self, wt):

        return wt*self.rhoL / self.Mm

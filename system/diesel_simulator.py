import numpy as np
from system.fluid import Fluid
from correlations.kinetics_diesel import rHDS_fun, rHDN_NB_fun, rHDN_B_fun, rHDA_fun, ft_reactants, effective_reactions
from correlations.bvp import OrthogonalCollocation
from scipy.integrate import solve_ivp
from scipy.optimize import root

class SolidConcentrations(Fluid):

    def __init__(self, cl0, heterogeneous = False, n_points = 3, ivp_tol=1e-16, method = "lm", **options):

        super().__init__(heterogeneous = heterogeneous, **options)

        self.cs0 = cl0
        self.method = method

        if self.heterogeneous:

            self.collocation = OrthogonalCollocation(
                    SolidConcentrations._transport_eq,
                    SolidConcentrations._bc_eq,
                    n_points, 3, x0=0, x1= self.reactor.dp/2)
            
        self.n_points = n_points
        self.ivp_tol = ivp_tol
        
        

    @staticmethod
    def _transport_eq(r, y, dy, d2y, yb, *args): #partial equations
        return d2y - ft_reactants(r, y, *args)
    
    @staticmethod    
    def _bc_eq(r, y, dy, d2y, yb, *args):  #boundary conditions
        return y - yb

    def mass_balance_surface(self, cs, cl, z):

        eff = np.ones(4)

        if self.heterogeneous:

            if z == 0:

                root_method = 'lm'
                cs[cs == 0] = self.ivp_tol
                y0 = np.column_stack((cs,) * (self.n_points + 1))

            else:
                y0 = self.collocation.y
                root_method = "hybr"

            args_ft = (cs, cl[7], self.T, self.rhoL, self.viscosity, self.vL, vH2, vH2S)
            self.collocation.collocate(y0, args=args_ft, method=root_method)

            args_reactions = (cl[7], self.T, self.rhoL)
            eff = self.collocation.effectiveness(effective_reactions, args_reactions)

        rHDS = rHDS_fun(cs,self.T)
        rHDN_NB = rHDN_NB_fun(cs, self.T, self.rhoL)
        rHDN_B = rHDN_B_fun(cs, self.T, self.rhoL)
        rHDA = rHDA_fun(cs, cl[7], self.T, self.rhoL)

        F = np.empty((7))

        F[0] =  1*eff[0] * self.rhob * rHDS  -  self.kSaS_oil*(cl[0] - cs[0])
        F[1] =  1*eff[1] * rHDN_NB -  self.kSaS_oil * (cl[1] - cs[1])
        F[2] =  -1*eff[2] * rHDN_B -  self.kSaS_oil * (cl[2] - cs[2])
        F[3] =  1*eff[3] * rHDA -  self.kSaS_oil * (cl[3] - cs[3])
        F[4] =  15*eff[0] * self.rhob * rHDS +6*eff[1] * rHDN_NB + 2*eff[2]*rHDN_B+3*eff[3] * rHDA   - self.kSaS_H2 * (cl[4] - cs[4])
        F[5] =  -9 * eff[0] * self.rhob * rHDS   - self.kSaS_H2S * (cl[5] - cs[5])
        F[6] =  -1 * eff[3] * rHDA -  self.kSaS_oil * (cl[6] - cs[6])

        return F
    
    def get_surface_concentrations(self, cl,z):

        cs = root(self.mass_balance_surface, self.cs0, (cl,z), method = self.method).x
        self.cs0 = cs

        return cs

class Simulator(Fluid):

    def __init__(self, y0, n_points = 3,n_points_integration=100, heterogeneous = False, **options):

        super().__init__(heterogeneous = heterogeneous, **options)

        y0[0:4] = np.vectorize(self.wt_to_molar)(y0[0:4])
        self.y0 = y0
        self.concentrations = SolidConcentrations(y0[:-2], heterogeneous, n_points)
        self.sol = None
        self.z = self.reactor.z
        self.n_points_integration = n_points_integration
        self.cL0_profile = None
        self.cL1_profile = None
        self.cL2_profile = None
        self.cL3_profile = None
        self.cL4_profile = None
        self.cL5_profile = None
        self.cL6_profile = None
        self.p4G_profile = None
        self.p5G_profile = None

    def differential_equations(self,cl, pG4, pG5, cs):

        return np.array(
        [self.mass_balance_liquid(cl[0], cs[0]),
        self.mass_balance_liquid(cl[1], cs[1]),
        self.mass_balance_liquid(cl[2], cs[2]),
        self.mass_balance_liquid(cl[3], cs[3]),     
        self.mass_balance_gas_liquid_phase_H2(pG4, cl[4], cs[4]),
        self.mass_balance_gas_liquid_phase_H2S(pG5, cl[5], cs[5]),
        self.mass_balance_liquid(cl[6], cs[6]),
        self.mass_balance_gas_phase_H2(pG4, cl[4]),
        self.mass_balance_gas_phase_H2S(pG5, cl[5])])
    
    def dy(self, z, variables):

        pG4, pG5 = variables[-2:]
        
        solid_concentrations = self.concentrations.get_surface_concentrations(variables[:-1], z)

        return self.differential_equations(variables[:-2], pG4, pG5, solid_concentrations)
    
    def solve(self):

        t_span = [0, self.z]
        t_eval = np.linspace(0, self.z, self.n_points_integration)
        self.sol = solve_ivp(self.dy, t_span=t_span, y0=self.y0, t_eval=t_eval, method = "RK45")
        (self.cL0_profile,  self.cL1_profile, self.cL2_profile, self.cL3_profile, self.cL4_profile, self.cL5_profile, self.cL6_profile,
        self.p4G_profile, self.p5G_profile) = self.sol.y

    
    def get_surface_concentrations_profiles(self):

        surface_concentrations_profiles = np.apply_along_axis(
            lambda row: self.concentrations.get_surface_concentrations(row, 0) 
                        if np.array_equal(row, self.sol.y.T[0]) 
                        else self.concentrations.get_surface_concentrations(row, None),
            axis=1,
            arr=self.sol.y.T
        )
        return surface_concentrations_profiles.T


    def get_outlet_fluid_conditions(self):

        return self.sol.y[:,-1]
    
    def get_evectiveness_profile(self, surface_concentrations_profiles):


        effectiveness_profiles = np.empty((0,4))

        eff = np.ones(4)

        for i in range(self.n_points_integration):

            if self.heterogeneous:

                cs = surface_concentrations_profiles.T[i]

                if i == 0:

                        root_method = 'lm'
                        cs[cs == 0] = self.concentrations.ivp_tol
                        y0 = np.column_stack((cs,) * (self.concentrations.n_points + 1))

                else:
                    y0 = self.concentrations.collocation.y
                    root_method = "hybr"

                args_ft = (cs, self.p5G_profile[i], self.T, self.rhoL, self.viscosity, self.vL, vH2, vH2S)
                self.concentrations.collocation.collocate(y0, args=args_ft, method=root_method)

                args_reactions = (self.p5G_profile[i], self.T, self.rhoL)
                eff = self.concentrations.collocation.effectiveness(effective_reactions, args_reactions)

            effectiveness_profiles = np.vstack((effectiveness_profiles , eff))

        return effectiveness_profiles.T
from system.fluid import Fluid
from correlations.bvp import OrthogonalCollocation
from correlations.kinetics_diesel import rHDS_fun, rHDN_NB_fun, rHDN_B_fun, rHDA_fun, ft_reactants, effective_reactions


class SolidSolidConcentrations(Fluid):

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
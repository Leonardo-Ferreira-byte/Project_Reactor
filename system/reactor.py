import numpy as np
from correlations.mass_transfer import aS_fun

class Reactor:

    def __init__(self, z = 31.54, dp = 0.254,  reactor_diameter = 2.54):

        self.z = z
        self.reactor_diameter = reactor_diameter
        self.Ac = np.pi * self.reactor_diameter**2 / 4
        self.V_reactor = self.Ac * z
        self.dp = dp
        self.porosity = self.get_bed_fraction(reactor_diameter, dp)
        self.aS = aS_fun(self.porosity, dp)

    @staticmethod
    def get_bed_fraction(reactor_diameter, dp):

        return 0.38 + 0.073*(1 + (reactor_diameter/dp - 2)**2/(reactor_diameter/dp)**2)
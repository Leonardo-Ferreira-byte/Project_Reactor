from reactor.solubility import get_viscosity
from reactor.mass_transfer import aS_fun, kiL_aL_fun, kiS_aS_fun
from reactor.solubility import get_diffusion_oil, get_diffusion_H2, get_diffusion_H2S
from reactor.data import roW, vH2, vH2S, get_API, get_volume_molar_oil


def test_diffusion_coefficients():

    T = 370 + 273.15
    T_MeABP = (451 + 273.15)*1.8
    Mm1 = 420
    rho0 = 0.9146
    d15_6 = rho0 / roW
    API=22
    P = 10
    T = 370 + 273.15
    GL = 0.00572
    API = 22
    porosity = 0.48 
    dp = 0.172
    viscosity = get_viscosity(T, API)
    vL = get_volume_molar_oil(T_MeABP, d15_6, Mm1)
    D1 = 3.25E-5
    D2 = 1.33e-4
    D4 = 1.10e-4
    aS = aS_fun(porosity, dp)

    assert abs(kiL_aL_fun(D2, viscosity, rho0, P, T, GL) - 7.06e-3) < 1e-6
    assert abs(kiL_aL_fun(D4, viscosity, rho0, P, T, GL) - 6.42e-3) < 1e-6

    assert abs(kiS_aS_fun(D1, viscosity, rho0, P, T, GL, aS) - 2.73e-2) < 1e-4
    assert abs(kiS_aS_fun(D2, viscosity, rho0, P, T, GL, aS) - 7.00e-2) < 1e-4
    assert abs(kiS_aS_fun(D4, viscosity, rho0, P, T, GL, aS) - 6.17e-2) < 1e-4
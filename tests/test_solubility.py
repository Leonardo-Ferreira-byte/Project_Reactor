from correlations.solubility import get_viscosity
from correlations.solubility import get_diffusion_oil, get_diffusion_H2, get_diffusion_H2S
from correlations.data import roW, vH2, vH2S, get_API, get_volume_molar_oil


def test_viscosity():

    #https://petrowiki.spe.org/Calculating_PVT_properties
    assert get_viscosity((120-32)/1.8+273.15, 37.9) - 0.0230 < 1e-4

def test_diffusion_coefficients():

    T = 370 + 273.15
    T_MeABP = (451 + 273.15)*1.8
    Mm1 = 420
    rho0 = 0.9146
    d15_6 = rho0 / roW
    API = get_API(d15_6)
    viscosity = get_viscosity(T, API)
    vL = get_volume_molar_oil(T_MeABP, d15_6, Mm1)
    assert abs(get_diffusion_oil(T, viscosity, vL) - 3.41e-5) < 1e-7
    assert abs(get_diffusion_H2(T, viscosity, vL, vH2) - 1.49e-4) < 1e-6
    assert abs(get_diffusion_H2S(T, viscosity, vL, vH2S) - 1.23e-4) < 1e-6
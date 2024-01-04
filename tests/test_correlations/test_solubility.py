from correlations.solubility import get_viscosity


def test_viscosity():

    # https://petrowiki.spe.org/Calculating_PVT_properties
    assert abs(get_viscosity((120-32)/1.8+273.15, 37.9) - 0.0230) < 1e-4


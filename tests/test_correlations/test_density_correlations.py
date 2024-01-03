from correlations.density_correlations import oil_density

def teste_oil_density():

    #Korsten and Hoffman, 1996
    assert abs(oil_density(0.9146, 0.101325, 323.15) - 0.8945) < 1e-4
    assert abs(oil_density(0.9146, 0.101325, 520/1.8) - 0.9146) <1e-4

    #Mederos and Anchieta, 2012
    assert abs(oil_density(0.873, 0.0781, 520/1.8) - 0.873) < 1e-3
    assert abs(oil_density(0.873, 0.0781, 520/1.8+5) - 0.869) < 1e-3

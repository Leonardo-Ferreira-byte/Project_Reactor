from correlations.data import get_API

def test_get_API():

    #Jarullah, Mujtaba, 2011
    assert abs(get_API(0.8558) - 33.84237) < 1e-6
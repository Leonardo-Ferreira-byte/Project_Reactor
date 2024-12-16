```
def get_viscosity(T, API):
    """
    Function to get the viscosity of oil (Glaso, 1980).

    Parameters
    ----------
    API: float or int
        API gravity of oil.

    T: int or float
        Temperature in K.

    Returns
    -------
    mi_L: float
        Viscosity of oil sulfite in g/(cm*s).
    """

    T = T * 1.8
    a = 10.313 * np.log10(T - 460) - 36.447
    return 3.141e10 * (T - 460) ** (-3.444) * ((np.log10(API)) ** a) / 100
```
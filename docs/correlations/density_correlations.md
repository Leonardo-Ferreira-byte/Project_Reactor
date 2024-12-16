
$$
\left(0.167 + \left(16.181 \cdot 10^{-0.0425 \cdot \rho_0}\right)\right) \cdot \frac{P}{1000} 
- 0.01 \cdot \left(0.299 + \left(263 \cdot 10^{-0.0603 \cdot \rho_0}\right)\right) \cdot \left(\frac{P}{1000}\right)^2
$$

$$ \sum_{i=0}^n i^2 = \frac{(n^2+n)(2n+1)}{6} $$


# Density Correlations Module

This module contains functions for calculating the corrections to oil density based on pressure and temperature, as well as the oil density itself. The calculations are based on standard conditions and use specific parameters for temperature and pressure adjustments.

# Functions
```
deltarhoP(rho0, P)

Description: Computes the correction to oil density based on pressure.

Parameters:

rho0: int or floatDensity at standard conditions (15.6°C, 101.3 kPa) in lb/ft³.

P: int or floatPressure in psia.

Returns:

deltarhoP: floatDensity correction with pressure in lb/ft³.
```

deltarhoT(rho0, P, T)

Description: Computes the correction to oil density based on temperature.

Parameters:

rho0: int or floatDensity at standard conditions (15.6°C, 101.3 kPa) in lb/ft³.

P: int or floatPressure in psia.

T: int or floatTemperature in °R.

Returns:

deltarhoT: floatDensity correction with temperature in lb/ft³.

oil_density(rho0, P, T)

Description: Computes the oil density as a function of pressure and temperature.

Parameters:

rho0: int or floatDensity at standard conditions (15.6°C, 101.3 kPa) in g/cm³.

P: int or floatPressure in Pa.

T: int or floatTemperature in K.

Returns:

oil_density: floatOil density corrected for temperature and pressure in g/cm³.

Notes

The rho0 input for oil_density is converted from g/cm³ to lb/ft³ internally.

Temperature T is converted from Kelvin to Rankine internally.

Pressure P is converted from Pascal to psia relative to atmospheric pressure using the imported Patm constant.

Dependencies

The module imports Patm from correlations.data for pressure conversions.
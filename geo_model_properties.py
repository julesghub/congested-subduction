# core UW bit
import UWGeodynamics as GEO
u = GEO.UnitRegistry

# +
gravity = 9.8 * u.meter / u.second**2
Tsurf   = 273.15 * u.degK
Tint    = 1573.0 * u.degK

kappa   = 1e-6   * u.meter**2 / u.second 

boxLength = 6000.0 * u.kilometer
boxHeight =  800.0 * u.kilometer
boxWidth  = 3000.0 * u.kilometer

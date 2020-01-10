# core UW bit
from UWGeodynamics.scaling import units as u

import geo_model_properties as modprop

alpha   = 3.0e-5 / u.degK

# Density is defined using the dimensionless relative density
mantle_density   = 3400 * u.kilogram / u.metre**3 # at surface
ref_density      = mantle_density * (modprop.Tint-modprop.Tsurf) * alpha

# Material properties
um = {
    'name'     : 'upper mantle',
    'index'    : 0,
    'viscosity': 1. * u.pascal * u.second * 1e20, 
    'density'  : 0.,
    'cohesion' : 1e3 * u.megapascal,
    'cohesion2': 1e3 * u.megapascal,
}

lm = {
    'name'     : 'lower mantle',
    'index'    : 1,
    'viscosity': 100. * u.pascal * u.second * 1e20, 
    'density'  : 0.,
    'cohesion' : 1e3 * u.megapascal,
    'cohesion2': 1e3 * u.megapascal,
}


# 80Ma oceanic lithosphere
subplate1 = {
    'name'     : 'oceanic plate 1',
    'index'    : 3,
    'viscosity':  1.000e+05 * u.pascal * u.second * 1e20, 
    'density'  : -3.388e-01,
    'cohesion' :  1.250e+01 * u.megapascal,
    'cohesion2':  6.250e+00 * u.megapascal,        
}

subplate1_phase = {
    'name'     : 'oceanic plate 1 after phase change',
    'index'    : 2,
    'viscosity':  1.000e+05 * u.pascal * u.second * 1e20, 
    'density'  :  1.107e+00,
    'cohesion' :  1.250e+01 * u.megapascal,
    'cohesion2':  6.250e+00 * u.megapascal,        
}

subplate2 = {
    'name'     : 'oceanic plate 2',
    'index'    : 4,
    'viscosity':  1.000e+05 * u.pascal * u.second * 1e20, 
    'density'  :  6.040e-01,
    'cohesion' :  6.744e+01  * u.megapascal,
    'cohesion2':  3.372e+01 * u.megapascal,      
}

          
subplate3 = {
    'name'     : 'oceanic plate 3',
    'index'    : 5,
    'viscosity': 1.930e+04 * u.pascal * u.second * 1e20, 
    'density'  : 3.849e-01,
    'cohesion' : 1.213e+02 * u.megapascal,
    'cohesion2': 6.066e+01 * u.megapascal,      
}

subplate4 = {
    'name'     : 'oceanic plate 4',
    'index'    : 6,
    'viscosity': 9.641e+01 * u.pascal * u.second * 1e20, 
    'density'  : 2.228e-01,
    'cohesion' : 1.000e+03 * u.megapascal,
    'cohesion2': 1.000e+03 * u.megapascal,      
}

# weak back arc material properties
backArc1 = {
    'name'     : 'backArc1',
    'index'    : 7,
    'viscosity':  4.978e+03 * u.pascal * u.second * 1e20,   
    'density'  : -1.198e+00,
    'cohesion' :  1.250e+01 * u.megapascal,   
    'cohesion2':  6.250e+00 * u.megapascal,
}
backArc2 = {
    'name'     : 'backArc2',
    'index'    : 8,
    'viscosity': 1.726e+02 * u.pascal * u.second * 1e20, 
    'density'  : 1.162e-01,
    'cohesion' : 2.500e+01 * u.megapascal,   
    'cohesion2': 1.250e+01 * u.megapascal,
}

# transitional lithosphere
trans1 = {
    'name'     : 'trans1',
    'index'    : 9,
    'viscosity':  5.000e+03 * u.pascal * u.second * 1e20, 
    'density'  : -1.977e+00,
    'cohesion' :  1.000e+03 * u.megapascal,
}
trans2 = {
    'name'     : 'trans2',
    'index'    : 10,
    'viscosity':  1.502e+03 * u.pascal * u.second * 1e20, 
    'density'  :  2.550e-01,
    'cohesion' :  1.000e+03 * u.megapascal,
}

# cratonic lithosphere
craton1 = {
    'name'     : 'craton1',
    'index'    : 11,
    'viscosity':  5.000e+03 * u.pascal * u.second * 1e20, 
    'density'  : -2.118e+00,
    'cohesion' :  1.000e+03 * u.megapascal,
}

craton2 = {
    'name'     : 'craton2',
    'index'    : 12,
    'viscosity':  1.675e+03 * u.pascal * u.second * 1e20, 
    'density'  :  2.533e-01,
    'cohesion' :  1.000e+03 * u.megapascal,
}

# assume ribbon and buoyant strip have cratonic material properties
ribbon = {
    'name'     : 'ribbon',
    'index'    : 13,
    'viscosity':  1e5 * u.pascal * u.second * 1e20, 
    'density'  : -2.11,
    'cohesion' : 1e3 * u.megapascal, 
}

buoyStrip = {
    'name'     : 'buoyStrip',
    'index'    : 14,
    'viscosity':  1e5 * u.pascal * u.second * 1e20,    # strong 
    'density'  : -2.11,                                # assume cratonic density
    'cohesion' : 1e3 * u.megapascal, 
}

# define material list
material_list = [ um, lm,
                 subplate1, subplate1_phase, subplate2, subplate3, subplate4, 
                 backArc1, backArc2, 
                 trans1, trans2, 
                 craton1, craton2, 
                 ribbon, 
                 buoyStrip ]

# core UW bit
import UWGeodynamics as GEO
u = GEO.UnitRegistry

import geo_model_properties as modprop

alpha   = 3.0e-5 / u.degK

# Density is defined using the dimensionless relative density
um_surf_density  = 3400 * u.kilogram / u.metre**3 # at surface
abs_density      = um_surf_density * (1 - alpha * (modprop.Tint - modprop.Tsurf))
ref_density      = um_surf_density * (modprop.Tint-modprop.Tsurf) * alpha

eclogite_surf_density = 3500 * u.kilogram / u.meter**3
eclogite_density = eclogite_surf_density * ( 1 - alpha * (modprop.Tint - modprop.Tsurf))
abs_density, ref_density, eclogite_density

# Material properties
um = {
    'name'     : 'Model', # really the 'upper mantle'
    'index'    : 0,
    'viscosity': 1. * u.pascal * u.second * 1e20, 
    'density'  : abs_density
}

lm = {
    'name'     : 'lower mantle',
    'index'    : 1,
    'viscosity': 100. * u.pascal * u.second * 1e20, 
    'density'  : abs_density
}


# 80Ma oceanic lithosphere
subplate1 = {
    'name'     : 'oceanic plate 1',
    'index'    : 3,
    'viscosity': 1.000e+25 * u.pascal * u.second, 
    'density'  : 3.22255e3 * u.kg / u.m**3,
    'cohesion' : 1.250e+01 * u.megapascal,
    'cohesion2': 6.250e+00 * u.megapascal,        
}

subplate1_phase = {
    'name'     : 'oceanic plate 1 after phase change',
    'index'    : 2,
    'viscosity': 1.000e+25 * u.pascal * u.second,
    'density'  : eclogite_density,
    'cohesion' : 1.250e+01 * u.megapascal,
    'cohesion2': 6.250e+00 * u.megapascal,        
}

subplate2 = {
    'name'     : 'oceanic plate 2',
    'index'    : 4,
    'viscosity':  1.000e+25 * u.pascal * u.second, 
    'density'  :  3.34755e3 * u.kg / u.m**3,
    'cohesion' :  6.744e+01 * u.megapascal,
    'cohesion2':  3.372e+01 * u.megapascal,      
}


subplate3 = {
    'name'     : 'oceanic plate 3',
    'index'    : 5,
    'viscosity': 1.92964e+24 * u.pascal * u.sec, 
    'density'  : 3.31851e3 * u.kg / u.m**3,
    'cohesion' : 1.21324e+02 * u.megapascal,
    'cohesion2': 6.06619e+01 * u.megapascal,      
}

subplate4 = {
    'name'     : 'oceanic plate 4',
    'index'    : 6,
    'viscosity': 9.64084e+21 * u.pascal * u.second,
    'density'  : 3.29701e+03 * u.kg / u.m**3,
    # no yielding    
}

# +
# cratonic lithosphere
craton1 = {
    'name'     : 'craton1',
    'index'    : 11,
    'viscosity': 5.00000e+23 * u.pascal * u.second,
    'density'  : 2.98662e+03  * u.kg / u.m**3,
    'cohesion' : 1.3000e+02 * u.megapascal,
    'cohesion2': 6.500e+01 * u.megapascal,
}

craton2 = {
    'name'     : 'craton2',
    'index'    : 12,
    'viscosity': 1.67484e+23 * u.pascal * u.second,
    'density'  : 3.30100e+03 * u.kilogram / u.m**3,
    'cohesion' : 2.91998e+02 * u.megapascal,
    'cohesion2': 1.45999e+02 * u.megapascal,}
# -

# transitional lithosphere
trans1 = {
    'name'     : 'trans1',
    'index'    : 9,
    'viscosity': 5.00000e+23 * u.pascal * u.second, 
    'density'  : 3.00529e+03 * u.kg / u.m**3,
    'cohesion' : 4.000e+01 * u.megapascal,
    'cohesion2': 2.000e+01 * u.megapascal,
}
trans2 = {
    'name'     : 'trans2',
    'index'    : 10,
    'viscosity': 1.50142e+23 * u.pascal * u.second,
    'density'  : 3.30122e+03 * u.kg / u.m**3,
    'cohesion' : 1.500e+02 * u.megapascal,
    'cohesion2': 7.500e+01 * u.megapascal,
}

# weak back arc material properties
backArc1 = {
    'name'     : 'backArc1',
    'index'    : 7,
    'viscosity': 4.97782e+23 * u.Pa * u.sec,   
    'density'  : 3.10862e+03 * u.kg/u.m**3,
    'cohesion' : 1.250e+01 * u.megapascal,   
    'cohesion2': 6.250e+00 * u.megapascal,
}
backArc2 = {
    'name'     : 'backArc2',
    'index'    : 8,
    'viscosity': 1.71934e+22 * u.pascal * u.second,
    'density'  : 3.28283e+03 * u.kg/u.m**3,
    'cohesion' : 2.500e+01 * u.megapascal,   
    'cohesion2': 1.250e+01 * u.megapascal,
}

# assume ribbon and buoyant strip have cratonic material properties
ribbon = {
    'name'     : 'ribbon',
    'index'    : 13,
    'viscosity': 1e25 * u.pascal * u.second,
    'density'  : 2.98662e+03  * u.kg / u.m**3,
}

buoyStrip = {
    'name'     : 'buoyStrip',
    'index'    : 14,
    'viscosity': 1e25 * u.pascal * u.second,    # strong 
    'density'  : 2.98662e+03  * u.kg / u.m**3,         # assume cratonic density
}

# define material list
material_list = [ um, lm,
                 subplate1, subplate1_phase, subplate2, subplate3, subplate4, 
                 backArc1, backArc2, 
                 trans1, trans2, 
                 craton1, craton2, 
                 ribbon, 
                 buoyStrip ]

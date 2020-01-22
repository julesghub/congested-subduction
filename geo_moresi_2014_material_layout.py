#!/usr/bin/env python
# coding: utf-8

# **Material Layout**

# In[1]:


# core UW bit
import UWGeodynamics as GEO
from UWGeodynamics import visualisation as vis
from UWGeodynamics.scaling import units as u
from UWGeodynamics.scaling import dimensionalise
from UWGeodynamics.scaling import non_dimensionalise as nd

# import underworld.functions
import underworld.function as fn

import math
import numpy as np
import os
import scipy


# In[3]:


import geo_model_properties as modprop
import relrho_geo_material_properties as matprop


# In[ ]:


exec(open("geo_moresi_2014_scaling.py").read())


# In[ ]:


# I assume here the origin is a the top, front, middle
# 'middle' being the slab hinge at top, front
pert = 0.3  # nondimensional pert 
slab_xStart = 2500. * u.kilometer
slab_dx = 3000.0 * u.kilometer  # was 7000 km in Moresi 2014
slab_dy =  100.0 * u.kilometer
slab_dz = 3000.0 * u.kilometer # this is the entire domain width
slab_layers = 4

slab_crust = 7.0 * u.kilometer


# In[ ]:


backarc_dx = 1200. * u.kilometer
backarc_dy =  100. * u.kilometer
backarc_xStart = slab_xStart - backarc_dx
backarc_layers = 2

trans_dx =  350. * u.kilometer
trans_dy =  100. * u.kilometer
trans_xStart = slab_xStart - backarc_dx - trans_dx
trans_layers = 2

craton_dx = 750. * u.kilometer
craton_dy = 150. * u.kilometer
craton_xStart = slab_xStart - backarc_dx - trans_dx - craton_dx
craton_layers = 2

ribbon_dx =  500. * u.kilometer
ribbon_dy =   50. * u.kilometer
ribbon_dz = 1500. * u.kilometer 
ribbon_xStart = slab_xStart + 500. * u.kilometer

bouyStrip_dx = 500. * u.kilometer
bouyStrip_dy =  50. * u.kilometer
bouyStrip_xStart = slab_xStart + slab_dx - bouyStrip_dx


# In[ ]:


# define coordinate uw.functions
fn_x = GEO.shapes.fn.input()[0]
fn_y = GEO.shapes.fn.input()[1]
fn_z = GEO.shapes.fn.input()[2]


# In[ ]:


def slabGeo(x, y, dx, dy,):
    slabShape = np.array([ (x,y), (x+dx,y), (x+dx,y-dy), (x,y-dy), (x-pert,y-dy-pert), (x-pert,y-pert) ])
    return GEO.shapes.Polygon(slabShape)


# In[ ]:


def backArcGeo(x, y, dx, dy,):
    backArcShape = np.array([ (x,y), (x+dx,y), (x+dx-pert/2.,y-dy), (x,y-dy)])
    return GEO.shapes.Polygon(backArcShape)


# In[ ]:


def boxGeo(x, y, dx, dy):
    boxShape = np.array([(x,y), (x+dx,y), (x+dx,y-dy), (x,y-dy)]) 
    return (GEO.shapes.Polygon(boxShape))


# In[ ]:


top50km = fn_y > nd(-50.*u.kilometer)
fn_top50kmGeo = fn.branching.conditional( ( (top50km,  matprop.subplate1['index']),
                                            (   True,  matprop.subplate2['index']),  
                                          )
                                        )

outside_ribbon = fn_z < nd(ribbon_dz)
fn_ribbonGeo = fn.branching.conditional( ( (outside_ribbon, fn_top50kmGeo  ),
                                           (          True, matprop.ribbon['index']), 
                                         )
                                       )  


# In[ ]:


# Set materials
conditions_2d = [ 
               ( fn_y < nd( -600.0 * 10**3 * u.meter ), matprop.lm['index']) ,

                # define ribbon gemetry - this is not 2D/3D auto safe
               ( boxGeo(nd(ribbon_xStart), 0., nd(ribbon_dx), nd(ribbon_dy)), matprop.ribbon['index']),

                # define buoyant strip gemetry
               ( boxGeo(nd(bouyStrip_xStart), 0., nd(bouyStrip_dx), nd(bouyStrip_dy)), matprop.buoyStrip['index']),

                # define slab geometery
               ( slabGeo(nd(slab_xStart), -0.*nd(slab_crust)         , nd(slab_dx), nd(slab_dy)/slab_layers), matprop.subplate1['index']),
               ( slabGeo(nd(slab_xStart), -1.*nd(slab_dy)/slab_layers, nd(slab_dx), nd(slab_dy)/slab_layers), matprop.subplate2['index']),
               ( slabGeo(nd(slab_xStart), -2.*nd(slab_dy)/slab_layers, nd(slab_dx), nd(slab_dy)/slab_layers), matprop.subplate3['index']),
               ( slabGeo(nd(slab_xStart), -3.*nd(slab_dy)/slab_layers, nd(slab_dx), nd(slab_dy)/slab_layers), matprop.subplate4['index']),
                # dedine back arc geometry
               ( backArcGeo(nd(backarc_xStart),  -0.*nd(backarc_dy)/backarc_layers, nd(backarc_dx), nd(backarc_dy)), matprop.backArc1['index']),
               ( backArcGeo(nd(backarc_xStart),  -1.*nd(backarc_dy)/backarc_layers, nd(backarc_dx), nd(backarc_dy)), matprop.backArc2['index']),
                # define transition gemetry
               ( boxGeo(nd(trans_xStart),  -0.*nd(trans_dy)/trans_layers, nd(trans_dx), nd(trans_dy)), matprop.trans1['index']),
               ( boxGeo(nd(trans_xStart),  -1.*nd(trans_dy)/trans_layers, nd(trans_dx), nd(trans_dy)), matprop.trans2['index']),
                # define craton gemetry
               ( boxGeo(nd(craton_xStart), -0.*nd(craton_dy)/craton_layers, nd(craton_dx), nd(craton_dy)), matprop.craton1['index']),
               ( boxGeo(nd(craton_xStart), -1.*nd(craton_dy)/craton_layers, nd(craton_dx), nd(craton_dy)), matprop.craton2['index']),
                # otherwise upper mantle!
                ( True, matprop.um['index']),
             ] 


# In[ ]:


conditions_3d = [ 
               ( fn_y < nd( -600.0 * 10**3 * u.meter ), matprop.lm['index']) ,

                # define ribbon gemetry - this is not 2D/3D auto safe
               ( boxGeo(nd(ribbon_xStart), 0., nd(ribbon_dx), nd(ribbon_dy)), fn_ribbonGeo),

                # define buoyant strip gemetry
               ( boxGeo(nd(bouyStrip_xStart), 0., nd(bouyStrip_dx), nd(bouyStrip_dy)), matprop.buoyStrip['index']),

                # define slab geometery
               ( slabGeo(nd(slab_xStart), -0.*nd(slab_crust)         , nd(slab_dx), nd(slab_dy)/slab_layers), matprop.subplate1['index']),
               ( slabGeo(nd(slab_xStart), -1.*nd(slab_dy)/slab_layers, nd(slab_dx), nd(slab_dy)/slab_layers), matprop.subplate2['index']),
               ( slabGeo(nd(slab_xStart), -2.*nd(slab_dy)/slab_layers, nd(slab_dx), nd(slab_dy)/slab_layers), matprop.subplate3['index']),
               ( slabGeo(nd(slab_xStart), -3.*nd(slab_dy)/slab_layers, nd(slab_dx), nd(slab_dy)/slab_layers), matprop.subplate4['index']),
                # dedine back arc geometry
               ( backArcGeo(nd(backarc_xStart),  -0.*nd(backarc_dy)/backarc_layers, nd(backarc_dx), nd(backarc_dy)), matprop.backArc1['index']),
               ( backArcGeo(nd(backarc_xStart),  -1.*nd(backarc_dy)/backarc_layers, nd(backarc_dx), nd(backarc_dy)), matprop.backArc2['index']),
                # define transition gemetry
               ( boxGeo(nd(trans_xStart),  -0.*nd(trans_dy)/trans_layers, nd(trans_dx), nd(trans_dy)), matprop.trans1['index']),
               ( boxGeo(nd(trans_xStart),  -1.*nd(trans_dy)/trans_layers, nd(trans_dx), nd(trans_dy)), matprop.trans2['index']),
                # define craton gemetry
               ( boxGeo(nd(craton_xStart), -0.*nd(craton_dy)/craton_layers, nd(craton_dx), nd(craton_dy)), matprop.craton1['index']),
               ( boxGeo(nd(craton_xStart), -1.*nd(craton_dy)/craton_layers, nd(craton_dx), nd(craton_dy)), matprop.craton2['index']),
                # otherwise upper mantle!
                ( True, matprop.um['index']),
             ] 


# In[ ]:





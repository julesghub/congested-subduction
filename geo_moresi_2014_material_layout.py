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
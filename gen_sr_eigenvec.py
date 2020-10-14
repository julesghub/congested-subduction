# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# #### This file is for post model analysis.
# The strain rate tensor of tracer swarm in the model is extraced from the output resuults and converted to it's eigenvector components which are then saved to hdf5.
#
#

# %%
import underworld as uw
from underworld import function as fn

# %%
import h5py

# %%
mesh = uw.mesh.FeMesh_Cartesian(elementRes=(128,48,48),
                                minCoord=(0,0,0),
                                maxCoord=(1,1,1))

# %%
mesh.load('./mesh.h5')

# %%
'''
IDEA 1
open every velocityField file.
from it create the cell-centered tensor field
save it and open with XDMF to play

IDEA 2
open every velocityField file.
from it create the cell-centered 3x principle strainrate tensor field
save it and open with XDMF to play

'''

# %%
v  = mesh.add_variable(nodeDofCount=mesh.dim)
sr = mesh.subMesh.add_variable(nodeDofCount=6)

e1 = mesh.subMesh.add_variable(nodeDofCount=3)
e2 = mesh.subMesh.add_variable(nodeDofCount=3)
e3 = mesh.subMesh.add_variable(nodeDofCount=3)

# %%
v.load('velocityField-13.h5')

# %%
# define some functions to get the strain rate tensor, symmetric part of velocity gradient tensor
fn_gradV = v.fn_gradient
fn_eps = fn.tensor.symmetric( fn_gradV )

# %%
gradV = fn_gradV.evaluate(mesh.subMesh)
eps = fn_eps.evaluate(mesh.subMesh)

# %%
gradV.shape, eps.shape

# %%
eps

# %%
sr.data[:] = fn_eps.evaluate(mesh.subMesh)

# %%
fH = sr.save('cell_sr.h5')
mH = mesh.save('mesh.h5')
sr.xdmf('cell_sr', fH, 'strinrate', mH, 'mesh')

# %%
import numpy as np

# %%
eps[0]


# %%
def vec_2_mat(v):
    return np.array([ [v[0], v[3], v[4]],
                      [v[3], v[1], v[5]],
                      [v[4], v[5], v[2]] ] )


# %%
eps_tensor = np.array(list(map(vec_2_mat, eps)))

# %%
eps_tensor.shape

# %%
w, v = np.linalg.eig(eps_tensor[:])

# w - the eigenvalues
# v - the eigenvectors

# %%
for i in range(v.shape[0]):
    out = w[i]*v[i]
    e1.data[i] = out[0]
    e2.data[i] = out[1]
    e3.data[i] = out[2]

# %%
out[0]

# %%
fH = e1.save('cell_e1')
e1.xdmf('e1', fH, 'e1', mH, 'mesh' )

fH = e2.save('cell_e2')
e2.xdmf('e2', fH, 'e2', mH, 'mesh' )

fH = e3.save('cell_e3')
e3.xdmf('e3', fH, 'e3', mH, 'mesh' )

# %%
e1.data[0]

# %%
e2.data[0]

# %%
### w, v, w*v

# %%
w[0] * v[:,0]

# %%
xx

# %%
w.shape, v.shape

# %%
np.broadcast_to(w,v.shape)

# %%
w[0]*v[0]

# %%
w[0]

# %%
v[0]

# %%
hm = v.T*w.T
hm.T

# %%
t = np.array([ [[1,1,1],[2,2,2],[3,3,3]], [[1,1,1],[2,2,2],[3,3,3]]])

# %%
l = np.array([[4,5,6], [4,5,6]])

# %%
l.shape,t.shape

# %%
l[0]*t[0]

# %%
l.T[:]

# %%
np.broadcast_arrays(l, [2])

# %%
x = np.array([[1], [2], [3]])
y = np.array([4, 5, 6])
b = np.broadcast_arrays(x, y)

# %%
x,y,b

# %%

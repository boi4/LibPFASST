import h5py
import sys
import getMesh
import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import pylab

filename     ='test.h5'
file         = h5py.File(filename,'r')

nr_fields    = (file.get('/input/problemdefinition/nr_fields'))[()][0][0]

Nx         = ( file.get('input/discretization/Nx')         )[()][0][0]
xleft      = ( file.get('input/problemdefinition/x_left')  )[()][0][0]
xright     = ( file.get('input/problemdefinition/x_right') )[()][0][0]

Ny         = ( file.get('input/discretization/Ny')        )[()][0][0]
yup        = ( file.get('input/problemdefinition/y_up')   )[()][0][0]
ydown      = ( file.get('input/problemdefinition/y_down') )[()][0][0]

XX, YY, dx, dy = getMesh.getMesh(Nx,Ny,xleft,xright,yup,ydown)

sol = file.get('input/problemdefinition/q_initial')
if sol==None:
    print "Could not find dataset 'solution' in selected HDF5 file. Now exiting."
    sys.exit()
sol = sol[...]
sol = sol.reshape(Ny,Nx,nr_fields)

fig = plt.figure()
ax  = fig.gca(projection='3d')
surf = ax.plot_surface(XX, YY, sol[:,:,0], cstride=1, rstride=1, cmap=cm.coolwarm, linewidth=0.0)
plt.show()
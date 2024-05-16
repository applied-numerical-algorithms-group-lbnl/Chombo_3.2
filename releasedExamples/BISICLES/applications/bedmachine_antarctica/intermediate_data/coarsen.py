"""

produce coarser resolution versions of 

bisicles_bedmachine_500m.nc

(1km, 2km, 4km, 8km)

for convenience

"""

from netCDF4 import Dataset
import numpy as np

def coarsen_1D(u):
    """
    coarsen u(1D) by a factor of two
    """

    n = len(u)
    uc = 0.5 * ( u[0:n-1:2] + u[1:n:2] )

    #print ('u[1] - u[0] = {} {}'.format(type(u[0]), u[1]-u[0]))
    return uc

def coarsen_2D(u):
    """
    coarsen u(2D) by a factor of two
    """

    n,m = u.shape

    uc = 0.25 * ( u[0:n-1:2, 0:m-1:2]
                + u[1:n:2,   0:m-1:2]
                + u[0:n-1:2, 1:m:2]
                + u[1:n:2,   1:m:2])

    return uc

def coarsenc(name, fine_name):


    v_names = ['thk','topg','umod','btrc','umodc']
    fine_nc = Dataset(fine_name,'r')
    coarse_nc = Dataset(name, 'w')
    
    x_fine = fine_nc.variables['x'][:]
    y_fine = fine_nc.variables['y'][:]
    x_coarse, y_coarse = coarsen_1D(x_fine), coarsen_1D(y_fine)

    xdim = coarse_nc.createDimension('x',size=len(x_coarse))
    ydim = coarse_nc.createDimension('y',size=len(y_coarse))
    
    xv = coarse_nc.createVariable('x','f8',('x'))
    xv[:] = x_coarse

    yv = coarse_nc.createVariable('y','f8',('y'))
    yv[:] = y_coarse
   
    
    for v in v_names:
        vv = coarse_nc.createVariable(v  ,'f8',('y','x'))
        vv[:,:] = coarsen_2D(fine_nc.variables[v][:,:])
    

    

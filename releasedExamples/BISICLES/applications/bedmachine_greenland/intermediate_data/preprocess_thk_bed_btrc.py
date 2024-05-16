#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:37:14 2019

produce initial data Greenland
from bedmachine and measures netcdf files.
Output is a netcdf file conatining the 
minimal BISICLES initial data: 

   x,y cartesian coordinates (EPSG 3413)
   topg (bedrock), 
   thk (ice thickness),
   umod (observed speed),
   umodc ( weight w(x,y) in the misfit f_m(x,y) =  w (|u_model| - |u_obs|)^2)
   btrc (initial guess for C in BISICLES inverse problem)
@author: stephen
"""

MASK_OCEAN=0
MASK_OPENLAND=1
MASK_GROUNDED=2
MASK_FLOATING=3

from netCDF4 import Dataset
import numpy as np
from  scipy.interpolate import RectBivariateSpline
from scipy import ndimage
from osgeo import ogr, osr
from coarsen import add_projection_attr_greenland
def cell_to_face_1D(x):
    n = len(x)
    xf = np.zeros(n+1)
    xf[1:n] = 0.5*(x[0:n-1]+x[1:n])
    xf[0] = x[0] - 0.5*(x[1]-x[0])
    xf[n] = x[n-1] + 0.5*(x[n-1] - x[n-2])
    return xf

def plot(x,y,z,zmin,zmax,name, step=4):
    import matplotlib.pyplot as plt
    print ('plotting ' + name)
    plt.figure(figsize=(6,6))
    plt.subplot(1,1,1,aspect='equal')
    plt.pcolormesh(cell_to_face_1D(x[::step]),cell_to_face_1D(y[::step]),
                   z[::step,::step],vmin=zmin,vmax=zmax,cmap='RdBu_r')
    plt.colorbar()
    plt.title(name)
    plt.savefig('{}.png'.format(name),dpi=600)

def zeros_2D(x,y):
    return np.zeros((len(y),len(x)))


def read_umod_mouginot(x,y,measures_nc):
    # compute |u| from the Mouginot data (450 m res)
    ncu = Dataset(measures_nc, 'r')
    xu = ncu.variables["x"][:]
    yu = np.flipud(ncu.variables["y"][:])
    vx = np.flipud( ncu.variables["VX"][:,:] ) # image data - top left to bottom right
    vy = np.flipud( ncu.variables["VY"][:,:] )
    umod = np.sqrt(vx**2 + vy**2)
    print ('max |u| = {} m/a '.format(np.max(umod)))
    print ('umod.shape = {} xu.shape = {} yu.shape = {}'.format(umod.shape,xu.shape,yu.shape))
    print ('interpolating ...')
    #interpolation of velocity data onto bedmachine grid
    spl = RectBivariateSpline(yu,xu,umod,kx=1,ky=1)    
    return spl(y,x) 


def patch_hole(x,y,lo,hi, thk, topg, usrf, mask):

    sx = np.logical_and(x > lo[0], x < hi[0])
    sy = np.logical_and(y < lo[1], y < hi[0])

    hole_n = np.sum( np.where(thk[sy,sx] > 10, 1 , 0) ) 
       

def patch_holes(x,y,thk, topg, usrf, mask):
    #holes are odd regions of thin ice
    print ('patching holes')
    topg_f = ndimage.minimum_filter(topg, 4)
    BIG_DIFF = 300.0 
    SMALL_THK = 10.0
    
    hole = np.logical_and(topg_f < topg - BIG_DIFF, mask == MASK_GROUNDED)
    hole = np.logical_and(hole, thk < SMALL_THK)
    
    print (' found {} hole cells'.format(np.sum(np.where(hole,1,0))))
        
    topg_f = np.where(hole, topg_f, topg)
    thk = np.where(hole, usrf-topg_f,thk) 

    plot(x,y,topg_f-topg,-BIG_DIFF, BIG_DIFF, 'patch_filler', step=2)
    return thk,topg_f  

def remove_islands(thk,mask):
    from skimage import morphology
    print ('removing islands')
    eps = 1.0
    ice_free = thk < eps
    small_area = 128*128
    ice_free = morphology.remove_small_holes(ice_free, small_area, in_place = True) 

    thk = np.where(ice_free, 0, thk)
    
    return thk

def preprocess(output_nc, bedmachine_nc, measures_nc):


    C_MAX = 1.0e+4 # maximum value for C
    C_MIN = 1.0e+1 # minimum value for C

    #desired dimensions
    nx = 5120*2
    ny = 9216*2
    
    ncbm = Dataset(bedmachine_nc, 'r')
    xbm = ncbm.variables["x"][:]
    ybm = np.flipud(ncbm.variables["y"][:])
    
    nxbm,nybm = len(xbm),len(ybm)
    topg = np.zeros((ny, nx))
    thk = np.zeros((ny, nx))
    usrf_bm = np.zeros((ny, nx))
    mask = np.zeros((ny, nx))
    umod = np.zeros((ny, nx))

    print ('xbm[0] = {}, ybm[0] = {}'.format(xbm[0],ybm[0]))

    #bed machine data dimensions
    dx = xbm[1] - xbm[0]
    
    #desired data dimensions
    tol = 1.0e-10
    x = np.arange(xbm[0],xbm[0]+nx*dx+tol,dx)
    y = np.arange(ybm[0],ybm[0]+ny*dx+tol,dx)


    #bedmachine data
    topg[0:nybm,0:nxbm] =  np.flipud(ncbm.variables["bed"][:,:])
    thk[0:nybm,0:nxbm]  =  np.flipud(ncbm.variables["thickness"][:,:])
    usrf_bm[0:nybm,0:nxbm]  =  np.flipud(ncbm.variables["surface"][:,:])
    mask[0:nybm,0:nxbm]  = np.flipud(ncbm.variables["mask"][:,:])

    #speed data
    umod[0:nybm,0:nxbm]  = read_umod_mouginot(xbm,ybm,measures_nc)

    #thickness/bedrock mods 
    thk = remove_islands(thk,mask)
    
    #thk,topg = patch_holes(x,y,thk, topg, usrf_bm, mask)

    #raise ValueError('enough for now')
    
    #dependents
    eps = 1.0e-6
    rhoi = 917.0
    rhoo = 1027.0
    sg = topg + thk
    sf = (1.0 - rhoi/rhoo)*thk
    
    grounded = np.logical_and( thk > eps, sg + eps > sf)
    usrf = np.where( grounded, sg, sf )
        
    
    print ('umod c ...')
    #umodc is the weight w(x,y) in the misfit f_m(x,y) =  w (|u_model| - |u_obs|)^2
    umodc = np.where(umod > 1.0, 1.0, 0.0)
    umodc = np.where(thk > 10.0, umodc, 0.0)

    #surface gradient
    print ('grad s ...')
    usrf = ndimage.gaussian_filter(usrf, 4) # smooth
    grads = zeros_2D(x,y)
    grads[1:ny-1,1:nx-1] = 0.5 / dx *  np.sqrt(
        (usrf[1:ny-1,0:nx-2] - usrf[1:ny-1,2:nx])**2 + 
        (usrf[0:ny-2,1:nx-1] - usrf[2:ny,1:nx-1])**2 )
 
    #initial guess for C
    print ('btrc...')
    btrc = rhoi * 9.81 * grads * thk / (umod + 1.0)
    btrc = np.where(umod > 1, btrc, C_MAX)
    btrc = np.where(btrc < C_MAX, btrc, C_MAX)
    btrc = np.where(btrc > C_MIN, btrc, C_MIN)
    #smooth with slippy bias
    print ('    ...filtering')
    btrcs = ndimage.minimum_filter(btrc, 8)
    btrcs = ndimage.gaussian_filter(btrcs, 32)
    btrc = np.where(btrc < btrcs, btrc, btrcs) # retain slippy spots
    
    #no ice value for C
    btrc = np.where(thk > 0, btrc, 100.0)
    
    #ouput netcdf
    print ('writing ...')
    ncout = Dataset(output_nc,'w')
    #dimensions
    xdim = ncout.createDimension('x',size=nx)
    ydim = ncout.createDimension('y',size=ny)
    #var defs

    xv = ncout.createVariable('x','f8',('x'))
    yv = ncout.createVariable('y','f8',('y'))

    add_projection_attr_greenland(ncout, xv, yv)

    def create2D(name):
        v = ncout.createVariable(name,'f8',('y','x'))
        v.setncattr('grid_mapping','crs')
        return v

    topgv = create2D('topg')
    thkv = create2D('thk')
    umodv = create2D('umod')
    umodcv = create2D('umodc')
    btrcv = create2D('btrc')
    
    #data
    xv[:] = x
    yv[:] = y
    topgv[:,:] = topg
    thkv[:,:] = thk
    umodv[:,:] = umod
    umodcv[:,:] = umodc
    btrcv[:,:] = btrc

    ncout.close()

    dx = x[1] - x[0]
    print( ' {} < x < {} '.format(np.min(x) - 0.5 * dx, np.max(x) + 0.5*dx))
    dy = y[1] - y[0]
    print( ' {} < y < {} '.format(np.min(y) - 0.5 * dy, np.max(y) + 0.5*dy))


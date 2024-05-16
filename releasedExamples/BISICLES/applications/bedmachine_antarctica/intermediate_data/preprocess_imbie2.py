#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 11:50:48 2020

@author: stephen
"""

import matplotlib.pyplot as plt

def imbie2_mask_nc(output_nc_file, imbie2_nc_file, geometry_nc_file):
    
    import numpy as np
    from netCDF4 import Dataset
    from  scipy.interpolate import RectBivariateSpline
    
    masknc = Dataset(imbie2_nc_file,'r')
    xm = masknc.variables['x'][:]
    ym = masknc.variables['y'][:]
    basin =  1.0* masknc.variables['basinNumber'][:,:]
   # plt.pcolormesh(xm,ym,basin)
    masknc.close()
    
    nc = Dataset(geometry_nc_file, 'r')
    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
    nc.close()
    
    spl = RectBivariateSpline(xm,ym,basin,kx=1,ky=1)
    basin_fine = spl(x,y)
    
    ncout = Dataset(output_nc_file,'w')
    xdim = ncout.createDimension('x',size=len(x))
    ydim = ncout.createDimension('y',size=len(y))
    xv = ncout.createVariable('x','f8',('x'))
    yv = ncout.createVariable('y','f8',('y'))
    xv[:] = x
    yv[:] = y
    
    bv = ncout.createVariable('basin','f8',('y','x'))
    bv[:,:] = basin_fine
    ncout.close()
    
    
    return None

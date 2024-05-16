#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 15:12:38 2022

@author: stephen
"""

from netCDF4 import Dataset
import numpy as np
from  scipy.interpolate import RectBivariateSpline
from scipy import ndimage
from osgeo import ogr, osr
from coarsen import add_projection_attr_greenland

def preprocess(output_nc):
    """
    
    interpolate the atrm, acab and ghf fields
    onto the 4800 mesh define by 'bedmachine_greenland_4800m.nc'
    
    Parameters
    ----------
    output_nc : string
        netcdf file name to write

    Returns
    -------
    None.

    """
    
    geo_file = 'greenland_bedmachine_4800m.nc'
    ncbm = Dataset(geo_file, 'r')
    x = ncbm.variables['x'][:]
    y = ncbm.variables['y'][:]

    def regrid(x, y, nc, var_name, x_name, y_name, x_off=0, y_off=0):
        xnc = nc.variables[x_name][:] + x_off
        ync = nc.variables[y_name][ :] + y_off
        try:
            znc = nc.variables[var_name][0,:,:]
        except:
            znc = nc.variables[var_name][:,:]
            
        znc = np.where(znc > 1e6, 0, znc)    # seesm to detect the no data (1e+36)
            
        print(xnc.shape, ync.shape, znc.shape)
        spl = RectBivariateSpline(ync,xnc,znc,kx=1,ky=1) 
        
            
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(1,2,1,aspect='equal')
        plt.pcolormesh(xnc,ync,znc)
        plt.subplot(1,2,2,aspect='equal')
        plt.pcolormesh(x,y,spl(y,x))
        
        
        
        
        return spl(y,x) 

    # from Searise Greenland_5km_v1.1.nc, http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
    stemp_nc = Dataset('../external_data/surftemp_searise_epsg3413.nc', 'r') 
    stemp = regrid(x,y, stemp_nc, 'Band1','x','y')
    stemp += 273.15 # conver to kelvin
    stemp = np.where(stemp > 273.15, 273.15, stemp)
    smb_nc = Dataset('../external_data/smb_searise_epsg3413.nc', 'r')# file from Steve Price
    acab = regrid(x,y, smb_nc, 'Band1','x','y')
        
    ghf_nc =  Dataset('../external_data/geothermal_heat_flow_map_10km_with_NGRIP.nc','r') # Coglan et al 
    ghf = regrid(x,y, ghf_nc, 'GHF','X','Y')
    #convert from m J s-1 m-2 to J a-1 m-2
    ghf *= 3600.*24.*365.*1.0e-3

    #ouput netcdf
    print ('writing ...')
    ncout = Dataset(output_nc,'w')
    #dimensions
    xdim = ncout.createDimension('x',size=len(x))
    ydim = ncout.createDimension('y',size=len(y))
    # var defs
    xv = ncout.createVariable('x','f8',('x'))  
    yv = ncout.createVariable('y','f8',('y'))

    add_projection_attr_greenland(ncout, xv, yv)
    def create2D(name):
       v = ncout.createVariable(name,'f8',('y','x'))
       v.setncattr('grid_mapping','crs')
       return v
 
    stempv = create2D('stemp')
    acabv = create2D('acab')
    ghfv = create2D('ghf')

    
    #data
    xv[:] = x
    yv[:] = y
    stempv[:,:] = stemp
    acabv[:,:] = acab
    ghfv[:,:] = ghf


    ncout.close()
    return 
    
#zz = preprocess('test.nc') 
    
    
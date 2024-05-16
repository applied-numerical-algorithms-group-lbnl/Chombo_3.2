#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 09:26:36 2019

read an ISMIP6 2D nx,ny,nt  file and produce nt nc file
and hdf5 files

@author: steph
"""
NCTOAMR = '/home/stephen/Development/BISICLES-ISMIP6-AIS/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.OPT.ex'

NCTOAMR = '/global/homes/c/cornford/cori-bisicles/BISICLES-ismip6ais/code/filetools/nctoamr2d.Linux.64.CC.ftn.DEBUG.OPT.ex'

def split_times(nc_file_name, var_name, out_file_name_base, time_offset):
    """
    
    Read an ISMIP6 2D nx,ny,nt netcdf file and write
    a sequence of layered (nx,ny) netcdf files ready
    for conversion to hdf5 with the BISICLES nctoamr tool
    
    Parameters
    ----------

    nc_file_name : string
        input file name
    var_name : string
        name of the 3D variable to be extracted, ouput variables will be var_name_0000, etc
    out_file_name_base: string 
        first part of the ouput file name. a time index, and .nc will be appended  
    time_offset : int
        append time_offset + i to the name of file i in the time eequence

    """
    from netCDF4 import Dataset
    import numpy as np
    from  scipy.interpolate import RectBivariateSpline
    import os
    
    nc = Dataset(nc_file_name, 'r')
    x = nc.variables['x'][:].data
    y = nc.variables['y'][:].data
    t = nc.variables['time'][:]

    
    #nc file is on a node centered 761x761 8km mesh with a point on the south pole, 
    #but BISICLES AIS 8km mesh is 786x768 cell centred.
    
    xc = np.arange(4.0e+3,6144.0e+3,8.0e+3) - 6144.0e+3*0.5
    yc = xc.copy()
    
    
    for it,tt in enumerate(t):
        
        w_file_name = '{}_{:04d}'.format(out_file_name_base,it + time_offset) 
        nc_w_file_name = w_file_name + '.nc'
        ncw = Dataset(nc_w_file_name, 'w')
        
        xdim = ncw.createDimension('x',size=len(xc))
        xv = ncw.createVariable('x','f8',('x'))
        xv[:] = xc
        
        ydim = ncw.createDimension('y',size=len(yc))
        yv = ncw.createVariable('y','f8',('y'))
        yv[:] = yc
        
        zdim = ncw.createDimension('z',size=len(z))
        zv = ncw.createVariable('z','f8',('z'))
        zv[:] = z
        
        f2D = nc.variables[var_name][it,:,:,:]

        #set outer forcing to zero - should be well beyond shelves
        f2D = np.where(np.isnan(f2D),0.0,f_layer)
        # interpolate to cc grid
        spl = RectBivariateSpline(x,y,f2D,kx=1,ky=1)
        f_c = spl(xc,yc)
        fv = ncw.createVariable(var_name,'f8',('y','x'))
        fv[:,:] = f_c
            
        ncw.close()
      
        s = ''.join([' %s']*len(var_names)) % tuple(var_names)
        hdf5_w_file_name =  w_file_name + '.2d.hdf5'
        os.system('{} {} {} {}'.format( NCTOAMR, nc_w_file_name, hdf5_w_file_name, s))


#default arguments...
nc_file_name = 'HadGEM2-ES_RCP85_thermal_forcing_8km_x_60m.nc'
var_name = 'thermal_forcing'
out_file_name_base = 'HadGEM2-ES_RCP85_thermal_forcing_8km'
time_offset = 1995

#arguments
import sys
print(sys.argv)

if len(sys.argv) == 5:
    nc_file_name = sys.argv[1]
    var_name = sys.argv[2]
    out_file_name_base = sys.argv[3]
    time_offset  = int(sys.argv[4])


split_times(nc_file_name, var_name, out_file_name_base, time_offset)    
#print (nc_file_name, var_name, out_file_name_base, time_offset)    





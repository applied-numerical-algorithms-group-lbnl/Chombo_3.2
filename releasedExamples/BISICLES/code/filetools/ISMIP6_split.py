#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 09:26:36 2019

read an ISMIP6 3D nx,ny,nt  file and produce nt 'poor-man multidim' nc file
and hdf5 file

@author: steph
"""
#NCTOAMR = '/home/stephen/Development/BISICLES-ISMIP6-AIS/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.OPT.ex'

#NCTOAMR = '/global/homes/c/cornford/cori-bisicles/BISICLES-ismip6ais/code/filetools/nctoamr2d.Linux.64.CC.ftn.DEBUG.OPT.ex'

# point to "installed" and maintained version on cori
NCTOAMR = '/global/common/software/m1041/BISICLES/haswell/bin/nctoamr2d.Linux.64.CC.ftn.OPT.MPI.PETSC.ex'

def split(nc_file_name, var_name, out_file_name_base, time_offset, scale):
    """
    
    Read an ISMIP6 3D nx,ny,nz,nt netcdf file and write
    a sequence of layered (nx,ny) netcdf files ready
    for conversion to hdf5 with the BISICLES nctoamr tool
    
    Parameters
    ----------

    nc_file_name : string
        input file name
    var_name : string
        name of the 3D variable to be 'layered', ouput variables will be var_name_0000, etc
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
    #    x = nc.variables['x'][:].data
    x = np.arange(-3040.0e+3,3040.0e+3 + 1.0,8.0e+3)
    #    y = nc.variables['y'][:].data
    y = x.copy()
    t = nc.variables['time'][:]

    if len(x) != 761:
        print (len(x), np.min(x), np.max(x))
        e = ValueError() 
        raise(e)
    
    #nc file is on a node centered 761x761 8km mesh with a point on the south pole, 
    #but BISICLES AIS 8km mesh is 786x768 cell centred.
    #with the bottom left hand corner
    # at -3072.e3
    
    #BEDMAP 2 mesh defintion; cell ecntred
    #ncols         6667
    #nrows         6667
    #xllcorner     -3333500
    #yllcorner     -3333500
    #cellsize      1000

    #BISICLES covers 263:6144-263+1, 263:6144-263+1 (FORTRAN index)
    #BISICLES covers 262:6144-262, 263:6144-262 (c index)
    #cell center x[262] = -3071000.0
    x_0  = -3071000.0 - 500.0
    xc = np.arange(x_0 + 4.0e+3, x_0 + 8.0e+3*768 ,8.0e+3) 
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
    
        f = nc.variables[var_name][it,:,:]
        var_names = list()
    
        
        #set nan forcing to zero - should be well beyond shelves
        f = np.where(np.isnan(f),0.0,f)
        # interpolate to cc grid
        spl = RectBivariateSpline(x,y,f,kx=1,ky=1)
        f_c = spl(xc,yc)
        f_name = var_name
        var_names.append(f_name)
        fv = ncw.createVariable(f_name,'f8',('y','x'))
        fv[:,:] = scale*f_c
            
        ncw.close()
      
        s = ''.join([' %s']*len(var_names)) % tuple(var_names)
        hdf5_w_file_name =  w_file_name + '.2d.hdf5'
        os.system('{} {} {} {}'.format( NCTOAMR, nc_w_file_name, hdf5_w_file_name, s))


#default arguments...
nc_file_name = ' '
var_name = ' '
out_file_name_base = ' '
time_offset = 1995

#arguments
import sys
print(sys.argv)

if len(sys.argv) == 5:
    nc_file_name = sys.argv[1]
    var_name = sys.argv[2]
    out_file_name_base = sys.argv[3]
    time_offset  = int(sys.argv[4])

scale = 1.0
if (var_name == 'smb_anomaly'):
    spyr = 365*24*3600
    scale = spyr/918.0 # ISMIP units in kg m^-2 s^-1
    
split(nc_file_name, var_name, out_file_name_base, time_offset, scale)    
#print (nc_file_name, var_name, out_file_name_base, time_offset)    





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 17:55:04 2020

comon functions and data for ASE 

@author: stephen
"""




def coarsen_1D(u):
    """
    coarsen u(1D) by a factor of two
    """

    n = len(u)
    uc = 0.5 * ( u[0:n-1:2] + u[1:n:2] )
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


def interp_bilinear(xf,yf,xc,yc,uc):
    """
    interpolate uc(yc,xc) onto the grid defined by xf,uf
    """
    from  scipy.interpolate import RectBivariateSpline
    return RectBivariateSpline(yc,xc,uc,kx=1,ky=1)(yf,xf)

    

def coarsenc(name, fine_name,  v_names = ['thk','topg','umod','btrc','umodc']):

    from netCDF4 import Dataset    
   
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
        
GDAL = False

def write_ais_nc(file_name, x , y, var_dict):
    
    from netCDF4 import Dataset

    if GDAL:
        from osgeo import ogr, osr
    #ouput netcdf
    print ('writing ' + file_name)
    ncout = Dataset(file_name,'w')

    #dimensions
    xdim = ncout.createDimension('x',size=len(x))
    ydim = ncout.createDimension('y',size=len(y))

    #var def
    crsv =  ncout.createVariable('crs','int')

    EPSG = 3031
    if GDAL:
        crs = osr.SpatialReference()
        crs.ImportFromEPSG(EPSG)
        crs_wkt = crs.ExportToWkt()
        ncout.setncattr('spatial_ref',crs_wkt)
        crsv.setncattr('crs_wkt',crs_wkt)
        
    ncout.setncattr('Conventions','CF-1.7') 
    crsv.setncattr('EPSG',int(EPSG))
    crsv.setncattr('grid_mapping_name','polar_stereographic')
    crsv.setncattr('latitude_of_projection_origin', -90.0)
    crsv.setncattr('straight_vertical_longitude_from_pole', 0.0)
    crsv.setncattr('scale_factor',1.0)
    crsv.setncattr('standard_parallel',-71.0)
    crsv.setncattr('false_easting',0.0)
    crsv.setncattr('false_northing',0.0)


    xv = ncout.createVariable('x','f8',('x'))
    xv.setncattr('standard_name','projection_x_coordinate')
    xv.setncattr('units','meter')
    
    yv = ncout.createVariable('y','f8',('y'))
    yv.setncattr('standard_name','projection_y_coordinate')
    yv.setncattr('units','meter')

    def create2D(name):
        v = ncout.createVariable(name,'f8',('y','x'))
        v.setncattr('grid_mapping','crs')
        return v

    for k in var_dict:
        
        vv = create2D(k)
        vv[:,:] = var_dict[k]

    xv[:] = x
    yv[:] = y

    ncout.close()


def subset_nc(x_range, y_range, var_names, file_name):
    
    from netCDF4 import Dataset 
    import numpy as np
    
    ncbm = Dataset(file_name, 'r')
    x = ncbm.variables['x'][:]
    y = ncbm.variables['y'][:]

    #find bounds
    x_lo,x_hi = x_range
    i = np.where(np.logical_and(x >= x_lo, x < x_hi))
    xs = slice(i[0][0],1+i[0][-1])
    y_lo,y_hi = y_range
    j = np.where(np.logical_and(y >= y_lo, y < y_hi))
    ys = slice(j[0][0],1+j[0][-1])

    #subset data needed
    x = x[xs]
    y = y[ys]
    
    z = dict()
    for var_name in var_names:
        z[var_name] = ncbm.variables[var_name][ys,xs]
    
    return x,y,z

def subset_bedmachine(x_range, y_range, bedmachine_file):
    from netCDF4 import Dataset 
    import numpy as np
    
    ncbm = Dataset(bedmachine_file, 'r')
    x = ncbm.variables['x'][:]
    y = ncbm.variables['y'][:]

    #find bounds
    x_lo,x_hi = x_range
    i = np.where(np.logical_and(x >= x_lo, x < x_hi))
    xs = slice(i[0][0],1+i[0][-1])
    y_lo,y_hi = y_range
    j = np.where(np.logical_and(y >= y_lo, y < y_hi))
    ys = slice(j[0][0],1+j[0][-1])

    #subset data needed
    x = x[xs]
    y = y[ys]
    umod = ncbm.variables['umod'][ys,xs]
    umodc = ncbm.variables['umodc'][ys,xs]
    btrc = ncbm.variables['btrc'][ys,xs]
    thk = ncbm.variables['thk'][ys,xs]
    topg = ncbm.variables['topg'][ys,xs]

    return x,y,thk,topg,umod,umodc,btrc
    


def subset_bedmachine_temperature(x_range, y_range, bedmachine_temperature_file, n_layer, prefix='temp',form='{}{:06d}'):
    from netCDF4 import Dataset 
    import numpy as np
    
    ncbm = Dataset(bedmachine_temperature_file, 'r')
    x = ncbm.variables['x'][:]
    y = ncbm.variables['y'][:]

    #find bounds
    x_lo,x_hi = x_range
    i = np.where(np.logical_and(x >= x_lo, x < x_hi))
    xs = slice(i[0][0],1+i[0][-1])
    y_lo,y_hi = y_range
    j = np.where(np.logical_and(y >= y_lo, y < y_hi))
    ys = slice(j[0][0],1+j[0][-1])

    #subset data needed
    x = x[xs]
    y = y[ys]
    
    Tdict = {} # pack the layers into a dictionary
    
    for layer in range(0,n_layer):
        name = form.format(prefix,layer)
        Tdict[name] = ncbm.variables[name][ys,xs]
    
    return x,y,Tdict



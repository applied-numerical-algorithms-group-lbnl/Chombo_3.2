#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 17:19:20 2017

@author: stephen
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
from osgeo import gdal, osr, ogr
from osgeo.gdalconst import *
import os
def open_raster(raster_path):
    """
    Opens a tiff as specified by the user

    Returns a gdal dataset and an array of the raster
    """

    driver = gdal.GetDriverByName('Gtiff')
    driver.Register()
    src = gdal.Open(raster_path, GA_ReadOnly)
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    
    
    data=src.ReadAsArray()
    print("Opened %s" %(raster_path))
    
    return ulx,uly,xres,yres,np.flipud(data[:,:])

def save_raster(path,x,y,data):
    """
    wriye x,y,data to a geotiff path
    """

    driver = gdal.GetDriverByName('Gtiff')
    driver.Register()
    dx = x[1]-x[0]
    dataset = driver.Create(path, len(x), len(y), 1, gdal.GDT_Float64)
    dataset.SetGeoTransform( ( np.min(x), dx, 0, np.max(y), 0, -dx) ) 

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3031)
    wkt_projection = srs.ExportToWkt()
    dataset.SetProjection(wkt_projection)
    np.ma.set_fill_value(data,-9999.0)
    #data = np.flipud(data.filled())
    data = np.flipud(data)
    dataset.GetRasterBand(1).WriteArray(data)
    dataset.GetRasterBand(1).SetNoDataValue(-9999.0)
    dataset.FlushCache()  # Write to disk.


    
#%%
def save_nc(x,y,z_dict,filename):
    nc = Dataset(filename,'w')
    xdim = nc.createDimension('x',size=len(x))
    ydim = nc.createDimension('y',size=len(y))

    xvar = nc.createVariable('x','f8',('x'))
    yvar = nc.createVariable('y','f8',('y'))
    for k in z_dict:
        var  = nc.createVariable(k,'f8',('y','x'))   
        var[:,:] = z_dict[k]
        print (filename, k, np.min(var[:,:]), np.max(var[:,:]))
                
    xvar[:] = x
    yvar[:] = y

    nc.close()
        


#%%
def open_hogg_mask(file):  
    
    ox,oy,dx,dy,mask = open_raster(file)
    #ox,oy = -1.72855e+06,-158450.0
    eps = 1.0
    Ny,Nx = mask.shape
    #dx = 100.0
    xa = np.arange(ox,ox+dx*Nx-eps,dx )
    ya = np.flipud(np.arange(oy,oy-dx*Ny+eps,-dx ))
    splm = RectBivariateSpline(ya,xa,mask,kx=1,ky=1)
    maskb = splm(y,x)
    maskb = np.where(maskb > 0.99, 1.0, 0.0)
    #umoda = np.ma.masked_array(umoda, np.logical_not(maska == 1))
    return x,y,1.0 - maskb



def makehdf5(year,extra_mask=None):
    
    file = '../hogg_data/masks/pinei_20{}.tif'.format(year+1) # need to calve * before* the year end.
    xa,ya,mask =   open_hogg_mask(file)
    bflux = np.where(mask < 0.5, -1.0e+5, 0.0) # large melt in calved regions
    
    xx,yy = np.meshgrid(x-x[0],y-y[0])
    bflux = np.where(xx < 200.0e+3, 0, bflux )
    bflux = np.where(xx > 250.0e+3, 0, bflux )
    
    ncf = 'hogg_pig_calving_{}.nc'.format(year)
    save_nc(x,y,{ 'bflux' : bflux}, ncf)
    nctoamr = '~/Development/BISICLES/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.OPT.ex'
    hdf5 = 'hogg_pig_calving_{}.2d.hdf5'.format(year)

    cmd = nctoamr + ' ' + ncf + ' ' + hdf5 + ' bflux'    
    os.system(cmd)  

    return mask

#%%
asegeo = Dataset('../../ase_bedmachine_500m.nc','r')
x = asegeo.variables['x'][:]
y = asegeo.variables['y'][:]
uo = asegeo.variables['umod'][:,:]
uc = asegeo.variables['umodc'][:,:]


makehdf5(14)
makehdf5(15)
makehdf5(16)
makehdf5(17)
makehdf5(18)
makehdf5(19)
#makehdf5(20)

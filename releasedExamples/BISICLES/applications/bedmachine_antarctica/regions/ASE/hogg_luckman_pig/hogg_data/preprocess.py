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
def open_shp(shp_path):
    

    driver = ogr.GetDriverByName("ESRI Shapefile")
    src = driver.Open(shp_path, 0)
    layer = src.GetLayer()

    #for feature in layer:
    #    print ( f.ExportToJson())
    layer.ResetReading()
    ls = layer[0]
    
 
    
    print ( ls.ExportToJson())


    
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
def open_hogg_mask_100m(file):  
    
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


def open_hogg_100m(file):  
    
    ox,oy,dx,dy,umod = open_raster(file)
    mask = np.where(umod > 10,1,0)   
    #ox,oy = -1.72855e+06,-158450.0
    eps = 1.0
    Ny,Nx = umod.shape
    #dx = 100.0
    xa = np.arange(ox,ox+dx*Nx-eps,dx )
    ya = np.flipud(np.arange(oy,oy-dx*Ny+eps,-dx ))
    spl = RectBivariateSpline(ya,xa,umod,kx=1,ky=1)
    splm = RectBivariateSpline(ya,xa,mask,kx=1,ky=1)
    umoda = spl(y,x) 
    maska = splm(y,x)
    maska = np.where(uc > 0, maska, 0)
    maska = np.where(maska < 1, 0, maska)
    umoda = np.where(maska < 1, 0, umoda)
    #umoda = np.ma.masked_array(umoda, np.logical_not(maska == 1))
    return x,y,umoda,maska

def open_hogg_year(year, extra_mask=None):
    import glob

    f = 'speeds/PIG_*_{}????_{}????.coffsN_mag_DuFil_yrF_masked.gc.tiff'
    g = f.format(year,year)
    print (g)
    f = glob.glob(g)

    print ('n_files = {} '.format(len(f)))
    
    xa,ya,umoda,maska = open_hogg_100m(f[0])
    
    denom = maska

    tol = 1.0e-8
    
    for ff in f[1:]:
        xb,yb,umodb,maskb = open_hogg_100m(ff)
        umoda = np.where(maskb < 1.0-tol, umoda, umoda + umodb)
        denom += maskb
    
    umoda = np.where( denom > tol, umoda/ (denom + tol) ,0)
    maska = np.where(denom > tol, 1, 0)
    
    
    mf = 'masks/pinei_20{}.tif'.format(year)
    xb,yb,maskb = open_hogg_mask_100m(mf)
    
    maska = np.where(maskb < tol, 0, maska)
    
    #maska = np.where(maska < 1, 0, maska)
    #umoda = np.where(maska < 1, 0, umoda)
    
    if type(extra_mask) != type(None):
        maska = np.where(extra_mask < tol, 0 , maska )
 
            
    return xa,ya,umoda,maska
    
#%%

def makehdf5(year,extra_mask=None):
    xa,ya,umoda,maska =   open_hogg_year(year,extra_mask)
    ncf = 'hogg_pig_{}.nc'.format(year)
    save_nc(x,y,{ 'uo' : umoda , 'uc' : maska}, ncf)

    nctoamr = '~/Development/BISICLES/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.OPT.ex'
    hdf5 = 'hogg_pig_{}.2d.hdf5'.format(year)

    cmd = nctoamr + ' ' + ncf + ' ' + hdf5 + ' uo uc'    
    os.system(cmd)  

    return maska     

#%%
asegeo = Dataset('ase_bedmachine_500m.nc','r')
x = asegeo.variables['x'][:]
y = asegeo.variables['y'][:]
uo = asegeo.variables['umod'][:,:]
uc = asegeo.variables['umodc'][:,:]
#makehdf5(19)
#makehdf5(18)
#makehdf5(20)
#makehdf5(16)
#makehdf5(15)
makehdf5(17)
#use measures for 2014
nc15 = Dataset('hogg_pig_19.nc','r') # best coverage
ucm =  nc15.variables['uc'][:,:]
ncf= 'hogg_pig_14.nc'
uo = np.where(ucm > 0.1, uo, 0.0)
save_nc(x,y,{ 'uo' : uo , 'uc' : ucm}, ncf)
hdf5 = 'hogg_pig_14.2d.hdf5'
nctoamr = '~/Development/BISICLES/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.OPT.ex'
cmd = nctoamr + ' ' + ncf + ' ' + hdf5 + ' uo uc'    
os.system(cmd)  
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 17:53:40 2020

@author: stephen
"""

import os
import sys
sys.path.append(os.getcwd() + '/../../python')

import ais_bedmachine as ais
from ais_bedmachine import write_ais_nc, subset_bedmachine

import numpy as np
from osgeo import gdal, osr
from osgeo.gdalconst import *

def open_raster(raster_path):
    """
    Opens a tiff as specified by the user

    Returns a gdal dataset and an array of the raster
    """

    driver = gdal.GetDriverByName('Gtiff')
    driver.Register()
    dataset = gdal.Open(raster_path, GA_ReadOnly)
    data=dataset.ReadAsArray()
    print("Opened %s" %(raster_path))

    v = np.flipud(data[:,:])
    #v = data[:,:]
    return v
                  

x_lo,y_lo = -1632e+3,-1404e+3
x_hi,y_hi = x_lo + 960.0e+3 , y_lo + 768.0e+3

bedmachine_file = '../../intermediate_data/antarctica_bedmachine_500m.nc'
geometry_file_name = 'getz_bedmachine_500m.nc'

geometry_hdf5_file_name = 'getz_bedmachine_geometry_500m.2d.hdf5'
inverse_hdf5_file_name = 'getz_bedmachine_inverse_500m.2d.hdf5'

x,y,thk,topg,umod,umodc,btrc = subset_bedmachine(
    (x_lo,x_hi), (y_lo,y_hi), bedmachine_file)


#gdal_rasterize -l ANT_Basins_IMBIE2_v1.6 -a basin_no -tr 500.0 500.0 -a_nodata 0.0 -te -1631750.0 -1403750.0 -351750.0 -635750.0 -ot Float32 -of GTiff /home/stephen/Development/BISICLES/applications/bedmachine_antarctica/external_data/Rignot_ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp getz_rignot.tif
#gdal_rasterize -l ANT_Basins_IMBIE2_v1.6 -a basin_no -tr 500.0 500.0 -a_nodata 0.0 -te -1631750.0 -1403750.0 -671750 -635750.0 -ot Float32 -of GTiff /home/stephen/Development/BISICLES/applications/bedmachine_antarctica/external_data/Rignot_ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp getz_rignot.tif

basin_map = open_raster('getz_rignot.tif') # this one is *only* for getz


getz_mask = np.where(thk > 0, 1 , 0)
getz_mask = np.where(basin_map > 2, 0, getz_mask )
xx,yy = np.meshgrid(x,y)
getz_mask = np.where(yy > -692e+3,0,getz_mask)

import matplotlib.pyplot as plt


getz_grounded_mask = np.where(basin_map == 2, 1, 0)
#plt.pcolormesh(x,y,getz_grounded_mask)
x_max_getz =  -785000
#1plt.axvline(x_max_getz)
y_max_getz =  -692e+3
#plt.axhline(y_max_getz)
bb_cond = np.logical_and(xx < x_max_getz, yy < y_max_getz)

getz_mask = np.where(np.logical_and(bb_cond, basin_map == 0), 1, 0)
getz_mask = np.where(np.logical_or(getz_mask ==1, getz_grounded_mask==1 ),1,0)
plt.pcolormesh(x,y,getz_mask)

umod = np.where(getz_mask == 1, umod, 0)
umodc = np.where(getz_mask == 1, umodc, 0)
btrc = np.where(getz_mask == 1, btrc, 1.0e+5)


hf = np.where(topg < 0, -1027.0/917.0 * topg, 0)
hab = thk - hf
thk = np.where(np.logical_and(getz_mask==0,hab < 0), 0, thk)

write_ais_nc(geometry_file_name, x, y,
            {'thk':thk,'topg':topg,'umod':umod,'umodc':umodc,'btrc':btrc, 'rignot_basin':getz_grounded_mask, 'smask':getz_mask})



import os
def system(cmd):
    print(cmd)
    os.system(cmd)

NCTOAMR='~/Development/BISICLES/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.OPT.ex'

import os

system('{} {} {} umod umodc btrc'.format(NCTOAMR,
                                            geometry_file_name,
                                            inverse_hdf5_file_name))

system('{} {} {} thk topg '.format(NCTOAMR,
                                      geometry_file_name,
                                      geometry_hdf5_file_name))

rignot_hdf5_file_name = 'getz_bedmachine_rignot_500m.2d.hdf5'
system('{} {} {} rignot_basin '.format(NCTOAMR,
                                      geometry_file_name,
                                      rignot_hdf5_file_name))






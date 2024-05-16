#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 17:53:40 2020

@author: stephen
"""

import os
import sys
import numpy as np
from scipy import ndimage
sys.path.append(os.getcwd() + '/../../python')

from ais_bedmachine import write_ais_nc,subset_bedmachine,subset_bedmachine_temperature, interp_bilinear,subset_nc

x_lo,y_lo = 1.8e6, -0.64e6
x_hi,y_hi = x_lo + 1024.0e+3, y_lo+512.0e+3

bedmachine_file = '../../intermediate_data/antarctica_bedmachine_500m.nc'
bedmachine_temperature_file = '../../intermediate_data/antarctica_bedmachine_temperature_morlighem_4km_24.nc'
mask_file = '../../intermediate_data/antarctica_bedmachine_imbie2_basins_4km.nc'
KNOX_BASIN_NO = 3
geometry_file_name = 'knox_bedmachine_500m.nc'
temperature_file_name = 'knox_bedmachine_temperature_4km_24.nc'


geometry_hdf5_file_name = 'knox_bedmachine_geometry_500m.2d.hdf5'
inverse_hdf5_file_name = 'knox_bedmachine_inverse_500m.2d.hdf5'
temperature_hdf5_file_name = 'knox_bedmachine_temperature_4km_24.2d.hdf5'
basin_hdf5_file_name = 'knox_bedmachine_basin_500m.2d.hdf5'

x,y,thk,topg,umod,umodc,btrc = subset_bedmachine((x_lo,x_hi), (y_lo,y_hi), bedmachine_file)

xb,yb,zb = subset_nc((x_lo,x_hi), (y_lo,y_hi), ['basin'], mask_file)
#zb = np.where(np.abs(zb['basin'] - KNOX_BASIN_NO) < 0.0001,zb['basin'],0)
#no need to impose basin restrictions?
zb = np.where(np.abs(zb['basin'] - KNOX_BASIN_NO) < 1e6,zb['basin'],0)
basin = interp_bilinear(x,y,xb,yb,zb)
basin = ndimage.maximum_filter(basin, 16)
basin = ndimage.gaussian_filter(basin, 16)

basin = np.where(basin > 0.5, 1.0, 0.0)
umodc = np.where(basin > 0.5, umodc, 0.0)
umodc = ndimage.minimum_filter(umodc, 2)# shrink to aovid fronts, etc

BTRC_MAX = np.max(btrc)
btrc_basin = np.where(basin > 0.5, btrc, BTRC_MAX)
btrc_basin = ndimage.gaussian_filter(btrc_basin, 16) # smooth the transition at the basin edge
btrc = np.where(basin > 0.5, btrc, btrc_basin)

umod = np.where(basin > 0.5, umod, 0.0)
hab = thk - np.where(topg <0, -1027.0/917.0*topg, 0)
thk = np.where(np.logical_and(basin < 0.5, hab < 0), 0, thk) # get rid of shelves outside the basin


write_ais_nc(geometry_file_name, x, y,
            {'thk':thk,'topg':topg,'umod':umod,'umodc':umodc,'btrc':btrc,'knox':basin})

n_layer = 24
x,y,Tdict = subset_bedmachine_temperature( (x_lo,x_hi), (y_lo,y_hi), bedmachine_temperature_file, n_layer)
write_ais_nc(temperature_file_name,x,y,Tdict)

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

system('{} {} {} knox '.format(NCTOAMR,
                                      geometry_file_name,
                                      basin_hdf5_file_name))

args = ''
for i in range(0,n_layer):
    args += 'temp{:06d} '.format(i)
    
system('{} {} {} {}'.format(NCTOAMR,temperature_file_name,
       temperature_hdf5_file_name, args))



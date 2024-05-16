#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 10:48:36 2020

@author: scornford
"""
import os
import sys
sys.path.append(os.getcwd() + '/../../python')

from ais_bedmachine import write_ais_nc
from bisiclesIO import BisiclesData,  nctoamr

import numpy as np
from netCDF4 import Dataset


dx='1km'
out_base = 'bedmachine_antarctica_post_inverse-500m'
ncf = '{}_{}.{}'.format(out_base,dx,'nc')
hdf = '{}_{}.{}'.format(out_base,dx,'2d.hdf5')

INTERMEDIATE_DATA_PATH='../../intermediate_data'

def geometry_file(dx):
    return '{}/antarctica_bedmachine_{}.nc'.format(INTERMEDIATE_DATA_PATH,dx)

ncgeo = Dataset(geometry_file(dx),'r')
x = ncgeo.variables['x'][:]
y = ncgeo.variables['y'][:] 
ncgeo.close()
def origin(xc):
    return xc[0]-0.5*(xc[1]-xc[0])

ctrl =  BisiclesData('inverse-500m/ctrl.ant_bmach.04lev.000060000002.2d.hdf5',
                     level=4,origin=(origin(x),origin(y)),
                     croplo=(0.0,0.0),plot_file=False)    

c_one  = ctrl.beta
c_third = ctrl.beta * (1.0 + ctrl.speed**(2.0/3.0))
uf = 300.0
c_third_jreg_300 = c_third * ( ctrl.speed/uf + 1.0 )**(1.0/3.0) 
mucoef = ctrl.mucoef

write_ais_nc(ncf, ctrl.x,ctrl.y, {'c_one':c_one, 'c_third':c_third,
                                  'c_third_jreg_300':c_third_jreg_300,  'mucoef':mucoef})
nctoamr(ncf,hdf,'c_one c_third  c_third_jreg_300 mucoef')

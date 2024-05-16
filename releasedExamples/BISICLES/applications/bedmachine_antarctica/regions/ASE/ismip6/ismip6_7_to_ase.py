#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues Oct 13 2020

@author: Dan
"""

import os
import sys
import numpy as np
from scipy import ndimage
sys.path.append(os.getcwd() + '/../../../python')
from netCDF4 import Dataset
#from amrfile import io as amrio
#amrio.freeAll()
import glob

NCTOAMR = '/home/dan/code/BISICLES-public/code/filetools/nctoamr2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex'


from ais_bedmachine import write_ais_nc,subset_bedmachine,subset_bedmachine_temperature, interp_bilinear,subset_nc,coarsen_1D

def get_xy(xy_file):
    from netCDF4 import Dataset 
    import numpy as np
    ncxy = Dataset(xy_file, 'r')
    x = ncxy.variables['x'][:]
    y = ncxy.variables['y'][:]

    return x,y

def interp_and_write(x_out,y_out,basic_path, var_path):
    fullpath = basic_path + var_path
    ase_prefix = 'ase_'
    newpath = basic_path + ase_prefix + var_path
    #    print newpath
    
    os.system('mkdir -p {}'.format(newpath))

    files = sorted(glob.glob(fullpath + '*.nc'))[:]    
    for k,file in enumerate(files):
        print file

        #set up file and shape it
        edited_file = file.replace(fullpath, ase_prefix)                
        final_file = newpath + edited_file
        ncout = Dataset(final_file, 'w')

        xdim = ncout.createDimension('x',size=len(x_out))
        xv = ncout.createVariable('x','f8',('x'))
        xv[:] = x_out
        
        ydim = ncout.createDimension('y',size=len(y_out))
        yv = ncout.createVariable('y','f8',('y'))
        yv[:] = y_out

        #open full-AIS file
        read_nc = Dataset(file, 'r')
        x_in = read_nc.variables['x'][:]
        y_in = read_nc.variables['y'][:]
        var_names = list()
        
        vars_in = read_nc.variables
        for var_name in vars_in.keys():
            # if there's a z, just preserve it
            if (var_name == 'z'):
                print 'adding z'
                z_in = read_nc.variables[var_name][:]
                zdim = ncout.createDimension('z',size=len(z_in))
                zv = ncout.createVariable('z','f8',('z'))
                zv[:] = z_in
            elif (var_name == 'gamma0'):
                print 'adding gamma0'
                gamma_in = read_nc.variables[var_name][:]
                print gamma_in
                gamma_out = ncout.createVariable(var_name,'f8',())
                gamma_out = gamma_in
                print gamma_out
            elif ((var_name != 'x') and (var_name != 'y')):
                #print var_name
                var_in = read_nc.variables[var_name][:]
            
                var_out = interp_bilinear(x_out,y_out,x_in, y_in, var_in)
            
                var_names.append(var_name)
                fv = ncout.createVariable(var_name,'f8',('y','x'))
                fv[:,:] = var_out
            
        ncout.close()

        s = ''.join([' %s']*len(var_names)) % tuple(var_names)        
        var_hdf5_name = final_file.replace('.nc','.2d.hdf5')
        print var_hdf5_name
    
        os.system('{} {} {} {}'.format( NCTOAMR, final_file, var_hdf5_name, s))
    


print 'starting'
x_lo,y_lo = -1.838e6, -0.880e+6
x_hi,y_hi = x_lo + 896.0e+3, y_lo + 1024.0e+3

# 500m resolution mesh
bedmachine_500m_mesh_file = 'ase_bedmachine_500m_xy.nc'
# 4km mesh
#bedmachine_4km_mesh_file = 'ase_bedmachine_4km_xy.nc'

xfine,yfine = get_xy(bedmachine_500m_mesh_file)
#xcrse,ycrse = get_xy(bedmachine_4km_mesh_file)

# now use coarsen functions to go from 500m->1km->2km->4km
#500m->1km
x_1km = coarsen_1D(xfine)
y_1km = coarsen_1D(yfine)
#1km->2km
x_2km = coarsen_1D(x_1km)
y_2km = coarsen_1D(y_1km)
#2km->4km
xcrse = coarsen_1D(x_2km)
ycrse = coarsen_1D(y_2km)

print 'after x,y'
print xcrse
print ycrse

print 'smb forcing'
path = '/home/dan/people/ISMIP6/'
atm_path = 'ismip6_atmos/noresm1-m_rcp2.6/'
#atm_var = 'smb_anomaly'

interp_and_write(xcrse,ycrse,path, atm_path)

print 'ocean forcing'
ocean_path = 'ismip6_ocean/noresm1-m_rcp2.6/'

interp_and_write(xcrse,ycrse,path, ocean_path)

print 'ocean parameterizations'
ocean_parameterization_path = 'ismip6_ocean/parameterizations/'

interp_and_write(xcrse, ycrse, path, ocean_parameterization_path)

print 'imbie mask'
imbie_mask_path = 'ismip6_ocean/imbie2/'

interp_and_write(xcrse, ycrse, path, imbie_mask_path)

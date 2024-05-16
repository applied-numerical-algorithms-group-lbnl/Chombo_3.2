#!/usr/bin/python
# -*- coding: utf-8 -*-

# aggregate.py

# Reads in multiple plotfiles for a single variable (can either be netCDF4 files or HDF5 files). HDF5 files can be multi-level
# or single-level files, but if multilevel, this code assumes that level 0 is 8km. Plotfiles are then aggregated 
# into a single netCDF4 file for each variable.

import numpy as np
from netCDF4 import Dataset
from amrfile import io as amrio
import glob
from sys import argv
import os
import fnmatch

amrio.freeAll()

# usage ./aggregate.py <experiment_name> <path to extracted variables>

script = argv[0]
experiment = argv[1]
path_to_extracted_variables = argv[2]

secondsperyear = 31556926.0
daysperyear = 365.24
starting_year = 2011
icesheet_institution_model = '_AIS_CPOM-LBL_BISICLES_'

def read_data(plotfile, variable):

    # if the plotfile is a netCDF4 file 
    if plotfile.endswith('.nc'):
        nc = Dataset(plotfile, 'r')

        x = nc.variables['x'][:].data
        y = nc.variables['y'][:].data
        var = nc.variables[variable][:][:].data

        nc.close()

    # if the plotifle is an HDF5 file
    if plotfile.endswith('.hdf5'):

	level = 0

        amrid = amrio.load(plotfile)
        lo, hi = amrio.queryDomainCorners(amrid, level)
        iord = 1

        x, y, var = amrio.readBox2D(amrid, level, lo, hi, variable, iord)

        amrio.free(amrid)
    
    return x, y, var

def regrid_8km(orig_grid):
    # This function regrids an 8km BISICLES grid directly to an 8km ISMIP6 grid
    lb = 4 #lower bound
    ub = 765 #upper bound

    h = 0.25 * ( orig_grid[lb:ub, lb:ub] + orig_grid[lb+1:ub+1, lb:ub] + orig_grid[lb:ub, lb+1:ub+1] + orig_grid[lb+1:ub+1, lb+1:ub+1] )

    return h	


def create_nc_field(x, y, t, arr, variable, expt):
    
    # These variables have different names for the ISMIP6 submission
    var_dic = {'thickness': 'lithk', 'Z_base': 'topg', 'Z_bottom': 'base', 'Z_surface': 'orog', 'xVel': 'xvelmean', 'yVel': 'yvelmean', 'iceFrac': 'sftgif'}

    # If the variable we're interested in is part of the above dic,
    # then swap the name of the variable
    if variable in var_dic:
        variable = var_dic[variable]

    filename = variable + icesheet_institution_model + expt + '.nc'

    # Open an empty nc Dataset to write
    nc = Dataset(filename, 'w')

    # Create the dimensions
    xdim = nc.createDimension('x', size=len(x))
    ydim = nc.createDimension('y', size=len(y))
    tdim = nc.createDimension('time', size=len(t))

    # Create each variable
    xvar = nc.createVariable('x', 'f4', ('x'))
    yvar = nc.createVariable('y', 'f4', ('y'))
    tvar = nc.createVariable('time', 'f4', ('time'))
    ncvar = nc.createVariable(variable, 'f4', ('time', 'y', 'x'))

    # Set the variable attributes for time (ISMIP6 metadata standard)
    tvar.units = 'days_since 1-1-' + str(starting_year)
    tvar.calendar = '365.24_day'
    tvar.axis = 'T'
    tvar.long_name = "time"
    tvar.standard_name = 'time'

    # Set the variables
    xvar[:] = x
    yvar[:] = y

    # different handling of time depending on whether it is a flux variable or not
    if variable in ['acabf', 'libmassbffl', 'dlithkdt', 'licalvf', 'ligroundf']:
        tvar[:] = (t + 0.5) * daysperyear
    else:
        tvar[:] = t * daysperyear

    ncvar[:, :, :] = arr

    nc.close()
    
def aggregate(path, variable_name, expt_name):

    if variable_name in ['acabf', 'libmassbffl', 'dlithkdt', 'licalvf']:
    	num_of_plotfiles = len(fnmatch.filter(os.listdir(path + '/' + variable_name + '/' + expt_name + '/CF'), 'plot*CF*.hdf5'))
    else:    
	num_of_plotfiles = len(fnmatch.filter(os.listdir(path + '/' + variable_name + '/' + expt_name), 'plot*.hdf5'))

    t = np.arange(num_of_plotfiles)
    dx = 8 # resolution
    n = 761 # number of cells
    L = 3040.0e+3 #length from origin
    x = np.arange(-L, L+1, dx*1.0e+3)
    y = np.arange(-L, L+1, dx*1.0e+3)

    if variable_name in ['acabf', 'libmassbffl', 'dlithkdt', 'licalvf']:
	variable_files = sorted(glob.glob(path + '/' + variable_name + '/' + expt_name + '/CF/plot*CF*.hdf5'))[0:num_of_plotfiles]
    else:
	variable_files = sorted(glob.glob(path + '/' + variable_name + '/' + expt_name + '/plot*.hdf5'))[0:num_of_plotfiles]

    bulk = np.zeros( (len(t), n, n) )

    for k, file in enumerate(variable_files):
	print(file)
        xb, yb, varb = read_data(file, variable_name)
        bulk[k, :, :] = regrid_8km(varb)
    
    create_nc_field(x, y, t, bulk, variable_name, expt_name)

variables = ['thickness', 'Z_base', 'Z_surface', 
             'Z_bottom' , 'xVel'  , 'yVel',
             'iceFrac'  , 'acabf' , 'libmassbffl',
             'licalvf'  , 'ligroundf', 'sftgrf',
             'sftflf'   , 'strbasemag', 'dlithkdt']

for v in variables:
    print(v)
    aggregate(path_to_extracted_variables, v, experiment)

    

#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# calculateFluxScalars.py
#
# Create three separate netCDF4 files for tendacabf (S, FL), 
# tendlibmassbffl (S, FL), and tendlicalvf (S, FL) from the CF BISICLES 
# HDF5 plotfiles. These scalar values are created by spatially
# integrating over the entire grid of a single plotfile and then
# all files are aggregated into a single netCDF4 file.
#
# The code also takes in the flattened ligroundf (S, FL) produced 
# from running glfacesNew and creates a netCDF4 file for tendliground(S, FL). 
#
# The output from running this function will end up in the current working directory

import numpy as np
from netCDF4 import Dataset
from amrfile import io as amrio
import os
import glob
from sys import argv
import fnmatch

amrio.freeAll()

# usage: ./calculateScalars.py <path to flattened ligroundf output> <path to multi-component CF files> <experiment>

script = argv[0]
path_to_ligroundf = argv[1]
path_to_CF = argv[2]
experiment = argv[3]

secondsperyear = 31556926.0
daysperyear = 365.24
starting_year = 2011
icesheet_institution_model = '_AIS_CPOM-LBL_BISICLES_'

def load_CF_variables(plotfile, level):
    
    # Load the AMR BISICLES plotfile
    amrid = amrio.load(plotfile)
    lo, hi = amrio.queryDomainCorners(amrid, level)
    iord = 1

    # Extract the field variables acabf, libmassbffl, and licalvf
    x, y, acabf = amrio.readBox2D(amrid, level, lo, hi, 'acabf', iord)
    x, y, libmassbffl = amrio.readBox2D(amrid, level, lo, hi, 'libmassbffl', iord)
    x, y, licalvf = amrio.readBox2D(amrid, level, lo, hi, 'licalvf', iord)

    amrio.free(amrid)
    return x, y, acabf, libmassbffl, licalvf

def load_glfacesNew_variable(plotfile, level):
    # Note: the ligroundf variable is calculated using the glfacesNew filetool.
    # However, glfacesNew outputs multi-level BISICLES plotfiles, so these files
    # must be flattend to single-level 8km HDF5 first before running this function.

    # Load the flattened ligroundf plotfile
    amrid = amrio.load(plotfile)
    lo, hi = amrio.queryDomainCorners(amrid, level)
    iord = 1

    # Extract the ligroundf variable
    x, y, ligroundf = amrio.readBox2D(amrid, level, lo, hi, 'ligroundf', iord)

    amrio.free(amrid)
    return x, y, ligroundf

def create_nc_stat(t, arr, scalar, expt):

    filename = scalar + icesheet_institution_model + expt + '.nc'

    # Open an empty nc Dataset to write
    nc = Dataset(filename, 'w')

    # Create the time dimension
    tdim = nc.createDimension('time', size=len(t))

    # Create each variable (1 coordinate variable (time) and 1 regular variable)
    tvar = nc.createVariable('time', 'f4', ('time'))
    ncvar = nc.createVariable(scalar, 'f4', ('time'))

    # Set the variable attributes for time (ISMIP6 metadata standard)
    tvar.units = 'days_since 1-1-' + str(starting_year)
    tvar.calendar = '365.24_day'
    tvar.axis = 'T'
    tvar.long_name = "time"
    tvar.standard_name = 'time'

    # Set the variables
    tvar[:] = (t + 0.5) * daysperyear # 182.62, 
    ncvar[:] = arr

    nc.close()

def integrate(array):
    # Spatial integration over the field
    total = np.sum(array)
    return total

def calculate_ligroundf(path, expt):
# path is the directory where the flattened output of glfacesNew lives
    
    # Count the number of plotfiles in the directory
    num_of_plotfiles = len(fnmatch.filter(os.listdir(path), 'plot*.hdf5'))

    # 89 timesteps (this is different from the other CF variables because ligroundf
    # is generated from the 89 regular BISICLES plotfiles)
    t = np.arange(num_of_plotfiles)
    n = num_of_plotfiles

    # The files are single-level so level is 0
    level = 0

    # Create empty 1D arrays of size n
    tendligroundf = np.zeros(n)

    # Get a list of the flattened ligroundf files
    ligroundf_files = sorted(glob.glob(path + '/plot*.hdf5'))[0:num_of_plotfiles]

    # For each ligroundf file, grab the field values, integrate them, then
    # insert the scalar value into the empty array at index k. Each
    # kth file represents a timestamp in the experiment
    for k, file in enumerate(ligroundf_files):
        x, y, ligroundf = load_glfacesNew_variable(file, level)

        total_ligroundf = integrate(ligroundf)

        tendligroundf[k] = total_ligroundf
    
    create_nc_stat(t, tendligroundf, 'tendligroundf', expt)


def calculate_scalars_from_CF(CF_path, expt):
    
    #if expt == 'ismip6_ctrl2':
	#num_of_plotfiles = 119
    #else:
	#num_of_plotfiles = 90

    # Count the number of CF plotfiles in the directory
    num_of_plotfiles = len(fnmatch.filter(os.listdir(CF_path), 'plot*CF*.hdf5'))

    # 90 timesteps
    t = np.arange(num_of_plotfiles)
    n = num_of_plotfiles
    
    # The files are single-level so level is 0
    level = 0

    # Create empty 1D arrays of size n (to hold the scalar values)
    tendacabf = np.zeros(n)
    tendlibmassbffl = np.zeros(n)
    tendlicalvf = np.zeros(n)

    # Get a list of CF files in the given path
    CF_files = sorted(glob.glob(CF_path + '/plot*CF*.hdf5'))[0:num_of_plotfiles]

    # For each CF file, grab the field values, integrate them, then
    # insert the scalar value into the empty array at index k. Each
    # kth file represents a timestep in the experiment. 
    for k, file in enumerate(CF_files):
        x, y, acabf, libmassbffl, licalvf = load_CF_variables(file, level)
        
        total_acabf = integrate(acabf)
        total_libmassbffl = integrate(libmassbffl)
        total_licalvf = integrate(licalvf)

        tendacabf[k] = total_acabf
        tendlibmassbffl[k] = total_libmassbffl
        tendlicalvf[k] = total_licalvf
    
    # Create a netCDF4 file containing the time-aggregated scalar
    # array for each variable
    create_nc_stat(t, tendacabf, 'tendacabf', expt)
    create_nc_stat(t, tendlibmassbffl, 'tendlibmassbffl', expt)
    create_nc_stat(t, tendlicalvf, 'tendlicalvf', expt)

calculate_ligroundf(path_to_ligroundf, experiment)
calculate_scalars_from_CF(path_to_CF, experiment)    



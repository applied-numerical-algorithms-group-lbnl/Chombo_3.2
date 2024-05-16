#!/usr/bin/python 
# -*- coding: utf-8 -*-

# calculateStatsfromAllOut.py

# Read in the "<experiment_name>.allOut" file that was generated from the bisicles_stats-cori.py python script, 
# and create time-aggregated netCDF4 files for each variable (lim (S, ST), limnsw (S, ST), iareag (S, ST), iareaf (S, ST))

# example of a single line in <experiment_name>.allOut

# time = 2.010000000000e+03 iceVolumeAll = 2.626972156250e+16  iceVolumeAbove = 2.258467038158e+16 groundedArea = 1.209425400000e+13  floatingArea = 1.537522000000e+12  totalArea = 1.363177600000e+13 groundedPlusOpenLandArea = 1.214492700000e+13  iceMassAll = 2.411560439438e+19 iceMassAbove = 2.073272741029e+19  Total Melt = -2.192247381463e+11.

import re
import numpy as np
from netCDF4 import Dataset
from sys import argv
import os
import fnmatch

# usage: ./calculateStatsfromAllOut.py <ismip6_home_path> <path to .allOut file (including file)> <experiment name>

script = argv[0]
ismip6_home_path = argv[1]
allOut_path = argv[2]
experiment = argv[3]

secondsperyear = 31556926.0
daysperyear = 365.24
starting_year = 2011
icesheet_institution_model = '_AIS_CPOM-LBL_BISICLES_'

def create_nc_stat(t, arr, scalar, expt):

    filename = scalar + icesheet_institution_model  + expt + '.nc'

    # Open an empty nc Dataset to write
    nc = Dataset(filename, 'w')

    # Create the time dimension
    tdim = nc.createDimension('time', size=len(t))

    # Create each variable (1 coordinate variable (time) and 1 regular variable)
    tvar = nc.createVariable('time', 'f4', ('time'))
    ncvar = nc.createVariable(scalar, 'f4', ('time'))

    # Set the variable attributes for time
    tvar.units = 'days_since 1-1-'+ str(starting_year)
    tvar.calendar = '365.24_day'
    tvar.axis = 'T'
    tvar.long_name = "time"
    tvar.standard_name = 'time'

    # Set the variables
    tvar[:] = (t - starting_year) * daysperyear # days since the starting year
    ncvar[:] = arr

    nc.close()

def grab_scalars(line):

    time = float(re.search(r'time = (.*?) ', line).group(1))
    iceMassAll = float(re.search(r'iceMassAll = (.*?)  iceMassAbove', line).group(1))
    iceMassAbove = float(re.search(r'iceMassAbove = (.*?)  Total', line).group(1))
    groundedArea = float(re.search(r'groundedArea = (.*?)  floatingArea', line).group(1))
    floatingArea = float(re.search(r'floatingArea = (.*?)  totalArea', line).group(1))

    return time, iceMassAll, iceMassAbove, groundedArea, floatingArea


def calculate_stats_from_allOut(path, filename, expt):

    #if expt == 'ismip6_ctrl2':
	#num_of_plotfiles = 120
    #else:
	#num_of_plotfiles = 91

    # The number stats files that were created will tell us the number of timesteps
    num_of_statsfiles = len(fnmatch.filter(os.listdir(path + '/' + expt), 'plot*.stats'))

    # Number of timesteps
    n = num_of_statsfiles

    # Create empty 1-d arrays of size n for each variable
    time = np.zeros(n)
    lim = np.zeros(n)
    limnsw = np.zeros(n)
    iareag = np.zeros(n)
    iareaf = np.zeros(n)
    
    # Open the allOut file then read each line
    allOut = open(filename, 'r')
    experiment_data = allOut.readlines()

    # For each line, grab the scalar data from the line and insert it
    # into the kth index in the arrays created earlier
    for k, line in enumerate(experiment_data):
        data_per_timestep = grab_scalars(line)
        time[k], lim[k], limnsw[k], iareag[k], iareaf[k] = data_per_timestep

    # Create the netcdf files containing the filled arrays
    create_nc_stat(time, lim, 'lim', expt)
    create_nc_stat(time, limnsw, 'limnsw', expt)
    create_nc_stat(time, iareag, 'iareag', expt)
    create_nc_stat(time, iareaf, 'iareaf', expt)
    
    allOut.close()

calculate_stats_from_allOut(ismip6_home_path, allOut_path, experiment)


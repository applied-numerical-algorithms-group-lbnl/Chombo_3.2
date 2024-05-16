#!/usr/bin/python
# -*- coding: utf-8 -*-

# calculateTemperature.py

# Read in the first multi-level HDF5 BISICLES plotfile from an ISMIP6 experiment
# along with the temperature field file and grab variables litemptop (F), litempbotgr (F), 
# and litempbotfl (F) and creates three separate single-level netCDF4 files for each variable

import numpy as np
from amrfile import io as amrio
from netCDF4 import Dataset
import glob
from sys import argv

amrio.freeAll()

# usage: ./calculateTemperature.py <path to ismip6 experiment data> <path to temperature field file (including file)> <experiment_name>

script = argv[0]
experiment_path = argv[1]
temp_field = argv[2]
experiment_name = argv[3]

secondsperyear = 31556926.0
daysperyear = 365.24
starting_year = 2011
icesheet_institution_model = '_AIS_CPOM-LBL_BISICLES_'

def loadbike(plotfile, var, level):

    amrid = amrio.load(plotfile)
    lo, hi = amrio.queryDomainCorners(amrid, level)
    iord = 1
    x, y, v = amrio.readBox2D(amrid, level, lo, hi, var, iord)

    amrio.free(amrid)
    return x, y, v

def create_nc_field(x, y, t, arr, var, expt):

    filename = var + icesheet_institution_model + expt + '.nc'

    # open an empty nc Dataset to write
    nc = Dataset(filename, 'w')

    # Create the dimensions (this must be done before creating the variables)
    xdim = nc.createDimension('x', size=len(x))
    ydim = nc.createDimension('y', size=len(y))
    tdim = nc.createDimension('time', size =len(t))

    # Create each variable (3 coordinate variables and 1 regular variable)
    xvar = nc.createVariable('x', 'f4', ('x'))
    yvar = nc.createVariable('y', 'f4', ('y'))
    tvar = nc.createVariable('time', 'f4', ('time'))
    ncvar = nc.createVariable(var, 'f4', ('time', 'y', 'x'))

    # Set the attribute for time
    tvar.units = 'days_since 1-1-' + str(starting_year)
    tvar.calendar = '365.24_day'
    tvar.axis = 'T'
    tvar.long_name = "time"
    tvar.standard_name = 'time'

    # Set the variables
    xvar[:] = x
    yvar[:] = y
    tvar[:] = t * daysperyear
    ncvar[:,:,:] = arr

    nc.close()

def calculate_temperature(path, temperature_field, expt):

    level = 0 #Interpolation level
    n = 768 #BISICLES plotfile is 768x768
    t = np.arange(0.0, 1.0, 1.0) 

    # Grab the filename of the initial plotfile from a sorted list of plotfiles in an ISMIP6 experiment 
    initial_plotfile = sorted(glob.glob(path + '/plot*.hdf5'))[0]

    # Grab the filename of the temperature field - must be an 8km resolution temperature field
    #temperature_field = '/global/cscratch1/sd/cashafer/ISMIP6/scripts/antarctica-temperature-8km.2d.hdf5'
    #temperature_field = '/scratch/users/cashafer/ismip6/datafiles/results/antarctica-temperature-8km.2d.hdf5'

    # Load thickness and Z_base from the initial plotfile of a BISICLES ISMIP6 experiment (the first plotfile out of 89)

    x, y, thck = loadbike(initial_plotfile, 'thickness', level)
    x, y, topg = loadbike(initial_plotfile, 'Z_base', level)

    # Load the surface layer and base layer from the 8km temperature field

    x, y, surf_temp = loadbike(temperature_field, 'temp000000', level)
    x, y, base_temp = loadbike(temperature_field, 'temp000009', level)

    # Set some variables and calculate definitions for "grounded" and "shelf"

    rhoo = 1028.0 # specific ocean density
    rhoi = 918.0 # specific ice density

    hf = np.where(topg < 0.0, -topg * rhoo/rhoi , 0.0)
    hab = thck - hf
    eps = 1.0e-3
    epsice = 10.0
    ice = (thck > epsice)
    
    grounded = np.logical_and(hab > eps, ice)
    shelf = np.logical_and(ice, np.logical_not(grounded))

    # Define two more arrays bgr_temp (grounded basal temp) and bfl_temp (floating basal temp)
    # using base_temp defined above

    bgr_temp = np.where(grounded, base_temp, 0.0)
    bfl_temp = np.where(shelf, base_temp, 0.0)

    # Create three 1 x n x n arrays (in this case n = 768) to hold the temperature data for each variable

    bulk_surf = np.zeros( (1, n, n) )
    bulk_bgr = np.zeros( (1, n, n) )
    bulk_bfl = np.zeros( (1, n, n) )

    # Fill the arrays with the data

    bulk_surf[0,:,:] = surf_temp
    bulk_bgr[0,:,:] = bgr_temp
    bulk_bfl[0,:,:] = bfl_temp

    create_nc_field(x, y, t, bulk_surf, 'litemptop', expt)
    create_nc_field(x, y, t, bulk_bgr, 'litempbotgr', expt)
    create_nc_field(x, y, t, bulk_bfl, 'litempbotfl', expt)

calculate_temperature(experiment_path, temp_field, experiment_name)

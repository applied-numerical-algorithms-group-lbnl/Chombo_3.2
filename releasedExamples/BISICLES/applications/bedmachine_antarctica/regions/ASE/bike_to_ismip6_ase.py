# -*- coding: utf-8 -*-
"""
Extract ISMIP6 data from a directory of plot.X.2d.hdf5 files

PROTECT ASE version

Created Aug 11 2022

@author: ggslc
"""
import os
import numpy as np
from netCDF4 import Dataset
from amrfile import io as amrio
amrio.freeAll()
import glob

SECONDS_PER_YEAR = 31556926
DAYS_PER_YEAR = 365
YEARS_PER_SECOND = 1.0/SECONDS_PER_YEAR
TIME_START = 2008 #must be an integer !
ORIGIN_X = -1838000-500
ORIGIN_Y = -880000-500
PPISMIP6 = '~/Development/BISICLES/code/filetools/ppismip62d.Linux.64.g++.gfortran.DEBUG.OPT.ex'
ICE_DENSITY = 917.0
OCEAN_DENSITY= 1027.0
GRAVITY = 9.81

def load_bike(plotfile, var, level):
    """
    read a single variable from a bisicles plot.X.2d.hdf5 file
    
    Parameters
    ----------
    plotfile (string) : file name of the plot.X.2d.hdf5
    var  (string )    : name of the variable 
    level (int)       : copy/coarsen/refine data to this level

    Returns
    -------
    x (1D ndarray): x-coordinates of the mesh cell centres
    y (1D ndarray): y-coordinates of the mesh cell centres 
    v (2D ndarray): cell cenetred data v(y,x)
    t (float): time of the data

    """
    amrid = amrio.load(plotfile)
    lo,hi = amrio.queryDomainCorners(amrid, level)
    iord = 1
    x,y,v = amrio.readBox2D(amrid, level, lo, hi, var, iord)
    t = amrio.queryTime(amrid) - TIME_START
    amrio.free(amrid)
    
    return x, y, t, v

def bike_to_ismip6_ant(vbike,dx):
    """
    interpolate data 'vbike' from  bisicles 6144 km x 6144 km cell-centred
    antarctic domain to the smaller 6080 km x 6080 km ismip 6 node centred
    domain.
    
    Parameters:
        vbike (numpy ndarray): data to be interpolated
        dx (float): mesh resolution
    
    """
    n = int(6080.0 / dx) + 1 #ismip 6 grid
    m = int(6144.0 / dx) # bike grid
    k = int((m-n+1) / 2)             
    l = int(m-k+1)
    k -= 1
    l -= 1
    #print (k,l)
    h =  0.25 * ( vbike[k:l,k:l] + vbike[k+1:l+1,k:l] + vbike[k:l,k+1:l+1] + vbike[k+1:l+1,k+1:l+1])
   
    return h

def save_nc_ismip6(x,y,t,arr,var,units, model,expt, path):
    """ 
    Save a single uniform mesh to netcdf, one variable, several times 
    
    Parameters: 
        x,y (1D ndarray) : cell-centre co-ordinates
        t (1D ndarray)   : times in years 
        var (string)     : name of variable in netcdf file
        units(string)    : nam eof the units for the netcdf attribute
        arr(3D ndarray)  : data, in t,y,x order
        model(string)    : indicates the model config
        expt(string)     : experiment name
        path(string)     : path to save files
    """ 
    
    filename =  path + '/' + var + '_' + model + '_' + expt + '.nc'    
    
    nc = Dataset(filename, 'w')
    xdim = nc.createDimension('x',size=len(x))
    ydim = nc.createDimension('y',size=len(y))
    tdim = nc.createDimension('time',size=len(t))
   
    crs = nc.createVariable('crs','int64')
    crs.EPSG = 3031
    crs.grid_mapping_name = "polar_stereographic" ;
    crs.latitude_of_projection_origin = -90. ;
    crs.straight_vertical_longitude_from_pole = 0. ;
    crs.scale_factor = 1. ;
    crs.standard_parallel = -71. ;
    crs.false_easting = 0. ;
    crs.false_northing = 0. ;
    
    tvar  = nc.createVariable('time','f4',('time')) 
    tvar.units = 'days since 1-1-{:04d}'.format(TIME_START)
    tvar.calendar ='common_year'
    tvar.standard_name = "time"
    tvar.axis = 'time'
    tvar.long_name= "time"
    tvar[:] = t
        
    xvar = nc.createVariable('x','f4',('x'))
    xvar.standard_name = "projection_x_coordinate"
    xvar.units = 'meter'
    xvar[:] = x + ORIGIN_X
    
    yvar = nc.createVariable('y','f4',('y'))
    yvar.standard_name = "projection_y_coordinate"
    yvar.units = 'meter'
    yvar[:] = y + ORIGIN_Y
    
    ncvar =  nc.createVariable(var,'f4',('time','y','x'))
    ncvar.grid_mapping = "crs"
    ncvar.units = units
    ncvar[:,:,:] = arr
    
    
    nc.close()
    
def load_nc_imsip6(filename,var):
    """
    reload the data written in save_nc_ismip6(
    """
    nc = Dataset(filename, 'r')
    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
    v = nc.variables[var][:,:]
    nc.close()
    return x,y,v    





def create_complex_ismip6_files(bike_files, std_name, bike_names, function, 
                                units, level, model, expt, path):
    """
    generate compound variables from a list of bike_names and a function Sure this could be nicer....
    
    function should extract the data it needs through keyword args: a dictionary with members like 
    bike_name: array
    
    """
    n_time = len(bike_files)
    time = np.zeros(n_time)
    bulk = None 
    var_dict = dict()
    for k, file in enumerate(bike_files):
        
        for j, bike_name in enumerate(bike_names):
            x, y, t, var_dict[bike_name] = load_bike(file, bike_name, level)
            
        if type(bulk) is type(None):
            bulk = np.zeros((n_time, len(y) , len(x) ))
        
        bulk[k,:,:] = function(**var_dict)
        time[k] = t * DAYS_PER_YEAR # bike time in years
        
    save_nc_ismip6(x, y, time, bulk, std_name, units,  model,expt, path)

def create_simple_ismip6_files(bike_files, std_name, bike_name, scale, 
                               units, level, model, expt, path):
    """
    simple conversions, allowing for a change in name and scale
    
    """
    n_time = len(bike_files)
    time = np.zeros(n_time)
    bulk = None 
    for k, file in enumerate(bike_files) :
        x, y, t, v = load_bike(file, bike_name, level)
        if type(bulk) is type(None):
            bulk = np.zeros((n_time, len(y) , len(x) ))
        #print(file, t)
   
        time[k] = t * DAYS_PER_YEAR # bike time in years
        bulk[k,:,:] = v * scale

    save_nc_ismip6(x, y, time, bulk, std_name, units,  model,expt, path)


def create_ismip6_fluxes_fractions(bike_files, level, inputs_file, model, expt, path):
    """
    Create the ISMIP6 fluxes (ligroundf, licalvf) and fractions (sftgif,sftgrf,sftlif)
    
    Side effect: uses the ppismip6 tool to create an extra file 
    X.pp.2d.hdf5 for each plot file X.2d.hdf5
    
    need to define the inputs file

    """
    
    nmap = {'ligroundf': ('ligroundf', ICE_DENSITY * YEARS_PER_SECOND, 'kg m-1 s-1'),
            'licalvf': ('licalvf', ICE_DENSITY * YEARS_PER_SECOND, 'kg m-1 s-1'),
            'sftgrf': ('sftgrf', 1, '1'),
            'sftflf': ('sftflf', 1, '1'),}
    
    pp_files = bike_files.copy()
    
    for k, file in enumerate(bike_files):
        pp_files[k] = file.replace('.2d.hdf5','.pp.2d.hdf5',)
        if not os.path.exists(pp_files[k]):
            os.system(PPISMIP6 + ' ' + file + ' ' + inputs_file + ' ' + pp_files[k])
    
    for std_name, bike_name in nmap.items():
        print (std_name, bike_name)
        create_simple_ismip6_files(pp_files, std_name, bike_name[0], 
                                   bike_name[1], bike_name[2], level, model , expt, path)
    print ('land ice fraction')
    def sftgif(**kwargs):
        return kwargs['sftgrf'] + kwargs['sftflf']
    create_complex_ismip6_files(pp_files, 'sftgif', ['sftgrf','sftflf'], sftgif, '1', 
                                                  level, model , expt, path)
    

def create_ismip6_files(path, expt, model, inputs_file, level=0):
    """
    create the full set of ISMIP6 files from a directory of plot*hdf5 files
    
    """ 
   
    #snapshot data
    files = sorted(glob.glob(path + '/plot.' + expt + '.??????.2d.hdf5'))
     
    print('fluxes and fractions')
    create_ismip6_fluxes_fractions(files, level, inputs_file, model, expt, path)
    
    #simple data
    nmap = {'lithk': ('thickness', 1, 'm'), 
            'orog': ('Z_surface', 1, 'm'), 
            'topg': ('Z_base', 1, 'm'),
            'uvelmean':('xVel', YEARS_PER_SECOND, 'ms-1'),
            'vvelmean':('yVel', YEARS_PER_SECOND, 'ms-1')}
        
    for std_name, bike_name in nmap.items():
        print (std_name, bike_name)
        create_simple_ismip6_files(files, std_name, bike_name[0], bike_name[1], bike_name[2], level, model , expt, path)
    
    #basal friction. 
    print('strbasemag')

    def Tb(**kwargs):
        u = kwargs['xVel']
        v = kwargs['yVel']
        C = kwargs['dragCoef']
        h = kwargs['thickness']
        return C*np.sqrt(u*u + v*v)*h/(h+1e-10)
    
    
    create_complex_ismip6_files(files, 'strbasemag', ['xVel','yVel','dragCoef','thickness'], Tb, 'Pa', level, model , expt, path)



model = 'ASE_BISICLES_500m'    
expt='ase_fwd_bike_therm_ens24_uj50_iso_true_lev3'
inputs_file = 'inputs.' + expt
path = './fwd_bike_therm_ens'
inputs_file = path + '/inputs.' + expt
create_ismip6_files(path,expt,model,inputs_file, level=2)




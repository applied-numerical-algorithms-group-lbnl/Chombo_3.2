# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 12:59:17 2016

lim - land ice mass kg
limnsw - land ice mass not displacing sea water (vaf, kg)
iareagr - grouned area  (m^2)
iareafl - ice shelf area (m^2)
tendacabf - tendency of land ice mass due to SMB        (flux arosss surface?) kg/s
tendlibmassbf - tendency of land ice mass due to BMB     (flux across basse?) kg/s
tendlicalvf - tendency of land ice mass due to calving  (flux across CF) kg/s
tendligroundf - tendency of grounded ice mass thickness (flux across GL) kg/s


@author: ggslc
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from amrfile import io as amrio
import os
import glob
import csv
amrio.freeAll()


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

def loadgeo(plotfile, level, mask):
    print (plotfile)
    amrid = amrio.load(plotfile)
    lo,hi = amrio.queryDomainCorners(amrid, level)
    iord = 0
    x,y,thck = amrio.readBox2D(amrid, level, lo, hi, "thickness", iord)
    x,y,topg = amrio.readBox2D(amrid, level, lo, hi, "Z_base", iord) 
    x,y,acabf = amrio.readBox2D(amrid, level, lo, hi, "surfaceThicknessSource", iord)
    x,y,bmb = amrio.readBox2D(amrid, level, lo, hi, "activeBasalThicknessSource", iord) 
    hf = np.where(topg < 0.0, -topg* OCEAN_DENSITY/ICE_DENSITY , 0.0)
    hab = thck - hf
    t = amrio.queryTime(amrid) - TIME_START
    amrio.free(amrid)
    
    amrid = amrio.load(mask[0])
    x,y,mask = amrio.readBox2D(amrid, int(mask[2]), lo, hi, mask[1], 0)
    amrio.free(amrid)
    
    thck *= mask
    acabf *= mask
    bmb *= mask
    hab *= mask
    
    #plt.imshow(thck)
    
    return x,y,thck,topg,hab,acabf,bmb,t
 
def loadflux(plotfile, level):
    print (plotfile)
    amrid = amrio.load(plotfile)
    lo,hi = amrio.queryDomainCorners(amrid, level)
    iord = 0
    x,y,glflux = amrio.readBox2D(amrid, level, lo, hi, "liground", iord) 
    x,y,cfflux = amrio.readBox2D(amrid, level, lo, hi, "licalvf", iord) 
    amrio.free(amrid)
    return x,y,glflux,cfflux


def scalars(thck,topg,hab,smb,bmb,dx):
    
    lim = np.sum(thck)* ICE_DENSITY*dx*dx # mass, kg 
    hab = np.where(hab > 0.0, hab, 0.0)
    limnsw = np.sum(hab) * ICE_DENSITY*dx*dx #mass, kg
    
    eps = 1.0e-3
    ice = (thck > eps)
    grounded = np.logical_and(hab > eps, ice)
    shelf = np.logical_and( thck > eps, np.logical_not(grounded))
    iareag = dx*dx * np.sum(np.where(grounded, 1.0, 0.0))
    iareaf = dx*dx * np.sum(np.where(shelf, 1.0, 0.0))
    tendacabf = dx*dx * np.sum(np.where(ice,smb,0.0)) /SECONDS_PER_YEAR * ICE_DENSITY
    
    tendbflx = dx*dx * np.sum(np.where(ice,bmb,0.0)) /SECONDS_PER_YEAR * ICE_DENSITY
    
    gim =  dx*dx * np.sum(np.where(grounded,thck,0.0)) * ICE_DENSITY
    gmb = dx*dx * np.sum(np.where(grounded,smb,0.0)) /SECONDS_PER_YEAR* ICE_DENSITY
    
    
    return lim,limnsw,iareag,iareaf,tendacabf,tendbflx, gim, gmb


 
def add_nc_one_stat(nc,arr,var,unit):
    """ save a single variable to netcdf, one variable, several times """ 

    ncvar =  nc.createVariable(var,'f4',('time'))
    ncvar.units = unit
    ncvar[:] = arr

 
def save_nc_stats(scalars,model,expt,path):
 
    
      t,lim,limnsw,iareag,iareaf,tendacabf,tendlibmassf,tendlicalvf,tendligroundf = scalars
    
      filename =  path + '/scalar_' + model + '_' + expt + '.nc'    
    
      nc = Dataset(filename, 'w')
      tdim = nc.createDimension('time',size=len(t))
      tvar  = nc.createVariable('time','f4',('time')) 
      tvar[:] = t * DAYS_PER_YEAR
      tvar.units = 'days since 1-1-{:04d}'.format(TIME_START)
      tvar.calendar ='common_year'
      tvar.standard_name = "time"
      tvar.axis = 'time'
      tvar.long_name= "time"
   
      print ('(len(t), len(lim) )', len(t), len(lim))
      
      add_nc_one_stat(nc,lim,"lim","kg")
      add_nc_one_stat(nc,limnsw,"limnsw","kg")
      add_nc_one_stat(nc,iareag,"iareagr","m^2")
      add_nc_one_stat(nc,iareaf,"iareafl","m^2")
      add_nc_one_stat(nc,tendacabf,"tendacabf","kg s-1")
      add_nc_one_stat(nc,tendlibmassf,"tendlibmassbf","kg s-1")
      add_nc_one_stat(nc,tendlicalvf,"tendlicalvf","kg s-1")
      add_nc_one_stat(nc,tendligroundf,"tendligroundf","kg s-1")
      
      nc.close()

def load_nc_stats(ncfile):
     
     nc = Dataset(ncfile, 'r')
     #print (nc.variables)
     t = nc.variables['time'][:]
     lim =  nc.variables['lim'][:]
     limnsw =  nc.variables['limnsw'][:]
     iareag =    nc.variables['iareagr'][:]
     iareaf =    nc.variables['iareafl'][:]
     tendacabf =    nc.variables['tendacabf'][:]
     tendlibmassf =    nc.variables['tendlibmassbfl'][:]
     tendlicalvf =    nc.variables['tendlicalvf'][:]
     tendligroundf =    nc.variables['tendligroundf'][:]
     
     return t,lim,limnsw,iareag,iareaf,tendacabf,tendlibmassf,tendlicalvf,tendligroundf
     
def load_series(plotfiles, level = 3, mask = ('a',0)):
     
    n = len(plotfiles)
    tyr = np.zeros(n)
    lim = np.zeros(n)
    limnsw = np.zeros(n)
    iareag = np.zeros(n)
    iareaf = np.zeros(n)
    tendacabf = np.zeros(n)
    gim= np.zeros(n)
    gmb= np.zeros(n)
    tendligroundf = np.zeros(n)
    tendlibmassf = np.zeros(n)
    tendlicalvf = np.zeros(n)

    
    for k,file in enumerate(plotfiles) :
        x,y,thck,topg,hab,smb,bmb,tyr[k]  = loadgeo(file, level, mask)
        dx = x[1] - x[0]
        s = scalars(thck,topg,hab,smb,bmb,dx)
        lim[k],limnsw[k],iareag[k],iareaf[k],tendacabf[k],tendlibmassf[k],gim[k],gmb[k] = s
     
    tendligroundf[0] = (gim[1] - gim[0])/(tyr[1]-tyr[0])  
    tendligroundf[1:n-1] = (gim[2:n] - gim[0:n-2])/(tyr[2:n]-tyr[0:n-2])
    tendligroundf[n-1] = (gim[n-1] - gim[n-2])/(tyr[n-1]-tyr[n-2])
    tendligroundf = -tendligroundf/SECONDS_PER_YEAR + gmb
    tendligroundf *= -1 #flux inward +ve!
    
    tendlicalvf[0] = (lim[1] - lim[0])/(tyr[1]-tyr[0])  
    tendlicalvf[1:n-1] = (lim[2:n] - lim[0:n-2])/(tyr[2:n]-tyr[0:n-2])
    tendlicalvf[n-1] = (lim[n-1] - lim[n-2])/(tyr[n-1]-tyr[n-2])
    
    tendlicalvf = -tendlicalvf /SECONDS_PER_YEAR + tendacabf + tendlibmassf
    tendlicalvf *= -1 #flux inward +ve !

    
    return tyr,lim,limnsw,iareag,iareaf,tendacabf,tendlibmassf,tendlicalvf,tendligroundf


def read_plot_file_list(file_name, day = None):
    """
        
    read lists of plot files and times in the space-delimited format
    
    file_path year.day
    
    optionally subset on day

    Parameters
    ----------
    file_name(string) : name of the file 
    day : restrict to entries year.day == x.day

    Returns
    ------- 
    time (list) : time in *decimal* years
    plot_file (list)   : list of *hdf5* files  
    """

    time = list()
    plot_file = list()

    def append(year, day, file):
        #global time, plot_file
        time.append(float(year) + round(float(day)/365,2))
        plot_file.append(file)

    with open(file_name) as txt:
        txt_reader = csv.reader(txt, delimiter=' ') 
        for row in txt_reader:
            cur_year, cur_day = row[1].split(".")
            hdf5 = row[0].replace('nc','2d.hdf5')
            
            if (type(day) is int) and day != int(cur_day):
                i = 0 # ignore
            else:
                append(cur_year, cur_day, hdf5)
                    
    return time, plot_file        
        

m3sle = 1.0e-9  / 360 * 1.0e-3 #sse level equivanlent, 1 ton of ice
kgsle = m3sle / OCEAN_DENSITY

reload = True
mask = ('ase_bedmachine_basin_500m.2d.hdf5','ase',0)
model = 'ASE_BISICLES_500m' 
expt = 'XYZ'
path = 'fwd_measures_ocean-mask_basal-flux_5.fwd_ens14.356_uj150_lev3'

txt = path + '.times.txt'
time, plot_file =  read_plot_file_list(txt, day = int(182))

scalar_data = load_series(plot_file, level=3, mask=mask)
save_nc_stats(scalar_data, model, expt, path)  
    


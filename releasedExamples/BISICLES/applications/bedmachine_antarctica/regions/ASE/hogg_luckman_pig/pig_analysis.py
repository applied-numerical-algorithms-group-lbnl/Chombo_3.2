#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 17:16:30 2017

@author: stephen
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from osgeo import gdal, osr, ogr
from amrfile import io as amrio
import scipy.stats
import os
amrio.freeAll()


#%%

def save_raster(path,x,y,data):
    """
    wriye x,y,data to a geotiff path
    """

    driver = gdal.GetDriverByName('Gtiff')
    driver.Register()
    dx = x[1]-x[0]
    dataset = driver.Create(path, len(x), len(y), 1, gdal.GDT_Float64)
    dataset.SetGeoTransform( ( np.min(x), dx, 0, np.max(y), 0, -dx) ) 
    #dataset.SetGeoTransform( ( np.min(x), dx, 0, np.min(y), 0, dx) ) 
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3031)
    wkt_projection = srs.ExportToWkt()
    dataset.SetProjection(wkt_projection)
    np.ma.set_fill_value(data,-9999.0)
    #data = np.flipud(data.filled())
    data = np.flipud(data)
    dataset.GetRasterBand(1).WriteArray(data)
    dataset.GetRasterBand(1).SetNoDataValue(-9999.0)
    dataset.FlushCache()  # Write to disk.
#%%
def readplot(ctrlfile):
    #colormap of $\left.\tau_{nt}\right|_b$ and melt-rate ($m_i$)


    #frachi = 0.6 # don't want the whole domain
    amrid = amrio.load(ctrlfile)
    level=3
    lod,hid = amrio.queryDomainCorners(amrid, level)
    
    lo = [0,0]
    hi = [0,0]
    
    lo[0] = int(0.15*hid[0])
    hi[0] = int(0.35*hid[0])
     
    lo[1] = int(0.52*hid[1])
    hi[1] = int(0.68*hid[1])
    
    iord = 1
    x,y,thk = amrio.readBox2D(amrid, level, lo, hi, "thickness", iord)
    x,y,xvel = amrio.readBox2D(amrid, level, lo, hi, "xVelb", iord)
    x,y,yvel = amrio.readBox2D(amrid, level, lo, hi, "yVelb", iord)
    x,y,xvelo = amrio.readBox2D(amrid, level, lo, hi, "xVelo", iord)
    x,y,yvelo = amrio.readBox2D(amrid, level, lo, hi, "yVelo", iord)
    x,y,velc = amrio.readBox2D(amrid, level, lo, hi, "velc", iord)
    x,y,topg = amrio.readBox2D(amrid, level, lo, hi, "Z_base", iord)
    x,y,usrf = amrio.readBox2D(amrid, level, lo, hi, "Z_surface", iord)
    x,y,mucoef = amrio.readBox2D(amrid, level, lo, hi, "muCoef", iord)
    x,y,C = amrio.readBox2D(amrid, level, lo, hi, "Cwshelf", iord)
    hf = np.where(topg < 0.0, -topg*1028.0/918.0, 0.0)
    hab = thk - hf
    amrio.free(amrid)
    
    #amrid = amrio.load('getz-sectors-1km.2d.hdf5')
    #level=0
    #lo,hi = amrio.queryDomainCorners(amrid, level)
    #iord = 0
   # hi[0] = int(frachi*hi[0])
    #x,y,smask = amrio.readBox2D(amrid, level, lo, hi, "smask", iord)
    smask = None# np.where(smask == 3, 1, 0)
    #amrio.free(amrid)
    
    smask = np.where(thk > 0, 1,0)
    xo = -1831500*1.0e-3
    yo = -903500*1.0e-3
    return xo + x*(1e-3),yo + y*(1e-3),thk,xvel,yvel,xvelo,yvelo,velc,hab,usrf,topg,mucoef,C,smask


Cmax = 2.0e+4

def smask(d):
    return d[13]

def umodel(d):
    u = np.sqrt(d[4]**2 + d[3]**2)    + 1.0e-10
    u = np.ma.masked_array(u, smask(d) < 1)
    return u    
    
def uobs(d):
    uo =   np.sqrt(d[6]**2 + d[5]**2) + 1.0e-10
    uo = np.ma.masked_array(uo, d[7] < 0.5)
    
    return uo

def Cmodel(d):
    return  np.ma.masked_array(d[12], d[12] > Cmax)


def muCoef(d):
    return  np.ma.masked_array(d[11], d[12] > Cmax)

def Tbmodel(d):
    return (Cmodel(d) * umodel(d))*1.0e-3 #kPa

def hphimodel(d):
     return  muCoef(d) * d[2]

def ctof(x):
    N = len(x)+1
    xf = np.zeros(N)
    xf[1:N] = x+0.5*(x[1] - x[0])
    xf[0] = x[0] - 0.5*(x[1] - x[0])
    return xf   

def xkmf(d):
    
    return ctof(d[0])

def ykmf(d):
    return ctof(d[1])



def xkmc(d):
    
    return d[0]

def ykmc(d):
    return d[1]

def hab(d):
    return d[8]

#%% 
def umis(d):
    uo =  uobs(d)
    um =  umodel(d)
    #return np.ravel(uo) 
    mis = d[7] *(um-uo)#/um
    mis = np.ma.masked_array(mis, uo < 0.005)
    mis = np.ma.masked_array(mis, um < 0.005)
    mis = np.ma.masked_array(mis,d[7] < 0.5)
    return mis

matplotlib.rcParams.update({'font.size': 11})
matplotlib.rcParams.update({'font.family': 'serif'})

slist = [(1,1,1,0.1), (1,0.5,0.5,1), (1,1,0.0,1),(0.5,1,0.5,1), (0,1,1,1) , (0,0,1,1) , (1,0,1,1), (1,0,0,1) ]

steph_lin_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('toss',slist)

#%%

def finish(cblabel,cbticks=[-1,0,1],cbx=0.0,cbextend='both'):
     plt.xlabel('x (km)')
     plt.yticks([])
     cax = fig.add_axes([cbx, 0.85, 0.6 , 0.5])
     cax.get_xaxis().set_visible(False)
     cax.get_yaxis().set_visible(False)
     cax.set_frame_on(False)
     cb=plt.colorbar(aspect=10,orientation='horizontal', label=cblabel,shrink=0.5,extend=cbextend)
     cb.set_ticks(cbticks)
     plt.subplots_adjust(wspace=0.01,hspace=0.01,left=0.05,right=0.95,top=0.95,bottom=0.15)

#%%




#%%


def plotset(fileA,yearA,fileB,yearB):

    ctrlA =  readplot(fileA)
    ctrlB =  readplot(fileB)
    umin = 000.0
    umax = 5000.0
    labels = np.array( [0,1000,2000,3000,4000,5000] )
    ticks = labels 
    

    
    fig = plt.figure(figsize=(14,7))
    plt.subplot(241,aspect='equal')
    pcu = plt.pcolormesh(xkmf(ctrlA),ykmf(ctrlA),(umodel(ctrlA)),cmap=steph_lin_cmap ,vmin=umin,vmax=umax)
    plt.title('model speed ({})'.format(yearA))
    plt.ylabel('y (km)')
    plt.xticks([])
    plt.yticks([])
      
    plt.subplot(242,aspect='equal')
    pce = plt.pcolormesh(xkmf(ctrlA),ykmf(ctrlA),(umis(ctrlA)),cmap='RdBu_r' ,vmin=-200,vmax=200)
    plt.title('misift ({})'.format(yearA))
    plt.ylabel('y (km)')
    plt.xticks([])
    plt.yticks([])
    
    plt.subplot(243,aspect='equal')
    pcu = plt.pcolormesh(xkmf(ctrlB),ykmf(ctrlB),(umodel(ctrlB)),cmap=steph_lin_cmap ,vmin=umin,vmax=umax)
    plt.title('model. speed ({})'.format(yearB))
    plt.xticks([])
    plt.yticks([])
    
    plt.subplot(244,aspect='equal')
    pce = plt.pcolormesh(xkmf(ctrlA),ykmf(ctrlB),(umis(ctrlB)),cmap='RdBu_r' ,vmin=-200,vmax=200)
    plt.title('misift ({})'.format(yearB))
    plt.xticks([])
    plt.yticks([])
    
    plt.subplot(245,aspect='equal')
    pcdu = plt.pcolormesh(xkmf(ctrlB),ykmf(ctrlB),(umodel(ctrlB)-umodel(ctrlA)),cmap='RdBu_r' ,vmin=-500,vmax=500)
    plt.title('model speed change({}-{})'.format(yearB,yearA))
    plt.xticks([])
    plt.yticks([])
    
    
    plt.subplot(246,aspect='equal')
    pcphi = plt.pcolormesh(xkmf(ctrlA),ykmf(ctrlA),np.log2((muCoef(ctrlB)/muCoef(ctrlA))),cmap='PiYG_r' ,vmin=-0.25,vmax=0.25)
    plt.title(r'stiffening ({}/{})'.format(yearB,yearA))
    plt.xticks([])
    plt.yticks([])
    
    plt.subplot(247,aspect='equal')
    pcTb = plt.pcolormesh(xkmf(ctrlB),ykmf(ctrlB),np.log2(Tbmodel(ctrlB)/Tbmodel(ctrlA)),cmap='PiYG_r' ,vmin=-0.05,vmax=0.05)
    plt.title(r'basal drag ({}/{})'.format(yearB,yearA))
    plt.xticks([])
    plt.yticks([])
    
    plt.subplot(248,aspect='equal')
    pcC = plt.pcolormesh(xkmf(ctrlB),ykmf(ctrlB),np.log2(Cmodel(ctrlB)/Cmodel(ctrlA)),cmap='PiYG_r' ,vmin=-0.25,vmax=0.25)
    plt.title(r'basal drag coef. ({}/{})'.format(yearB,yearA))
    plt.xticks([])
    plt.yticks([])
    
    
    #finish(r'$\log \phi$',cbticks=[-1,0,1],cbx=0.0)
    
    plt.subplots_adjust(wspace=0.05,hspace=0.1,right=1.0,bottom=0.05,left=0.00,top=0.95)
    
    cbar_ax = fig.add_axes([0.02, 0.89 , 0.15 , 0.025])
    cb=plt.colorbar(pcu, cax=cbar_ax, orientation='horizontal',extend='max', ticks = ticks)     
    cb.set_label(r'Speed, $|u| $ (m/a)')
    
    cbar_ax = fig.add_axes([0.28, 0.89 , 0.15 , 0.025])
    cb=plt.colorbar(pce, cax=cbar_ax, orientation='horizontal',extend='both')      
    cb.set_label(r'$u_2 - u_1$ (m/a)')
    
    cbar_ax = fig.add_axes([0.79, 0.89 , 0.15 , 0.025])
    cb=plt.colorbar(pce, cax=cbar_ax, orientation='horizontal',extend='both')      
    cb.set_label(r'$u_2 - u_1$ (m/a)')
    
    cbar_ax = fig.add_axes([0.53, 0.89 , 0.15 , 0.025])
    cb=plt.colorbar(pcu, cax=cbar_ax, orientation='horizontal',extend='max', ticks = ticks)        
    cb.set_label(r'Speed, $|u| $ (m/a)')

    cbar_ax = fig.add_axes([0.02, 0.43 , 0.15 , 0.025])
    cb=plt.colorbar(pcdu, cax=cbar_ax, orientation='horizontal',extend='both')      
    cb.set_label(r'$u_2 - u_1$ (m/a)')

    cbar_ax = fig.add_axes([0.28, 0.43 , 0.15 , 0.025])
    cb=plt.colorbar(pcphi, cax=cbar_ax, orientation='horizontal',extend='both')      
    cb.set_label(r'$log_2(\phi_2/\phi_1)$')
       
    cbar_ax = fig.add_axes([0.53, 0.43 , 0.15 , 0.025])
    cb=plt.colorbar(pcTb, cax=cbar_ax, orientation='horizontal',extend='both')      
    cb.set_label(r'$log_2(\tau^b_2/\tau^b_1)$')
    
    cbar_ax = fig.add_axes([0.79, 0.43 , 0.15 , 0.025])
    cb=plt.colorbar(pcC, cax=cbar_ax, orientation='horizontal',extend='both')      
    cb.set_label(r'$log_2(\beta^2_2/\beta^2_1)$')

#%%
    
def plot_speeds(files,years):
    
    umin = 0.0
    umax = 0.5e+4
    nrow = len(files)
    index = np.arange(0,nrow)
    ncol = 6
    w_panel = 4.5
    h_panel = 3.5
    
    fig = plt.figure(figsize=(ncol*w_panel,nrow*h_panel))    
    j = 1
    
    u_prev = None
    Tb_prev = None
    D_prev = None
    
    for i, file,year in zip(index,files,years):
        ctrl =  readplot(file)
       
        u = umodel(ctrl)
        Tb =  Tbmodel(ctrl)
        D =  1.0 - muCoef(ctrl)
        du = umis(ctrl)
        plt.subplot(nrow, ncol, j ,aspect='equal') 
        j += 1
        pcu = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),u ,cmap=steph_lin_cmap ,vmin=umin,vmax=umax)
        plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
        plt.title('model speed ({})'.format(year))
        plt.xticks([])
        plt.yticks([])
        plt.colorbar(pcu)
        path='speed-{}.tif'.format(year)
        save_raster(path,xkmf(ctrl)*1.0e3,ykmf(ctrl)*1.0e3,u)
        
        plt.subplot(nrow, ncol, j ,aspect='equal') 
        j += 1
        pcu = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),uobs(ctrl) ,cmap=steph_lin_cmap ,vmin=umin,vmax=umax)
        plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
        plt.title('obs speed ({})'.format(year))
        plt.xticks([])
        plt.yticks([])
        plt.colorbar(pcu)
        
        plt.subplot(nrow, ncol, j ,aspect='equal') 
        j += 1
        pcu = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),du ,cmap='RdYlBu_r' ,vmin=-500,vmax=500)
        plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
        plt.title('misfit ({})'.format(year))
        plt.xticks([])
        plt.yticks([])
        plt.colorbar(pcu)
        
        

        if (i > 0):
            print('model speed up ({}-{})'.format(years[i-1],year))
            plt.subplot(nrow, ncol, j ,aspect='equal') 
            j += 1
            pcdu = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),u-u_prev ,cmap='RdYlBu_r' ,vmin=-500,vmax=500)
            plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
            plt.title('model speed up ({}-{})'.format(years[i-1],year))
            plt.xticks([])
            plt.yticks([])
            plt.colorbar(pcdu)
            
            #plt.subplot(nrow, ncol, j ,aspect='equal') 
            #j += 1
            #pcdTb = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),Tb-Tb_prev ,cmap='RdYlBu_r' ,vmin=-50,vmax=50)
            #plt.contour(xkmc(ctrl),ykmc(ctrl),Tbmodel(ctrl),[0.0e3])
            #plt.title(r'model $\Delta \tau_b$ (kPa) ({}-{})'.format(years[0],year))
            #plt.xticks([])
            #plt.yticks([])
            #plt.colorbar(pcdTb)
               
          
            
            plt.subplot(nrow, ncol, j ,aspect='equal') 
            j += 1
            dD = np.where(u > 1, D-D_prev, np.NaN)        
            pcdphi = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),dD ,cmap='Reds', vmin=0.0,vmax=0.5)
            path='delta_D-{}-{}.tif'.format(years[0],year)
            save_raster(path,xkmf(ctrl)*1.0e3,ykmf(ctrl)*1.0e3,np.where(dD>0,dD,0.0))
            plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
            plt.title(r'model $\Delta D$ ({}-{})'.format(years[0],year))
            plt.xticks([])
            plt.yticks([])
            plt.colorbar(pcdphi)
            
            #plt.subplot(nrow, ncol, j ,aspect='equal') 
            #j += 1
            #dTb = 2.0*(Tb-Tb_prev)/(Tb+Tb_prev)            
            #pcdtb = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),dTb ,cmap='RdYlBu_r', vmin=-0.5,vmax=0.5)
            #path='delta_Tb-{}-{}.tif'.format(years[0],year)
            #save_raster(path,xkmf(ctrl)*1.0e3,ykmf(ctrl)*1.0e3,dTb)
            #plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
            #plt.title(r'$2| \tau - \tau_0 | / | \tau + \tau _0| $ ({}-{})'.format(years[0],year))
            #plt.xticks([])
            #plt.yticks([])
            #plt.colorbar(pcdtb)
            
            plt.subplot(nrow, ncol, j ,aspect='equal') 
            j += 1         
            pcdphi = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),D ,cmap='Reds', vmin=0.0,vmax=1.0)
            path='D_{}.tif'.format(year)
            save_raster(path,xkmf(ctrl)*1.0e3,ykmf(ctrl)*1.0e3,np.where(D>0,D,0.0))
            
            plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
            plt.title(r'model $D$ ({})'.format(year))
            plt.xticks([])
            plt.yticks([])
            plt.colorbar(pcdphi)
            
        else:    
            
            #skip the speed up plot
            j += 1
            
            
            plt.subplot(nrow, ncol, j ,aspect='equal') 
            j += 1         
            pcdphi = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl),D ,cmap='Reds', vmin=0.0,vmax=1.0)
            path='D_{}.tif'.format(year)
            save_raster(path,xkmf(ctrl)*1.0e3,ykmf(ctrl)*1.0e3,np.where(D>0,D,0.0))
            
            plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
            plt.title(r'model $D$ ({})'.format(year))
            plt.xticks([])
            plt.yticks([])
            plt.colorbar(pcdphi)
            
            
            plt.subplot(nrow, ncol, j ,aspect='equal') 
            j += 1         
            pcTb = plt.pcolormesh(xkmf(ctrl),ykmf(ctrl), Tb ,cmap='RdYlBu_r', vmin=0.0,vmax=3.0e2)
            path='Tb_{}.tif'.format(year)
            save_raster(path,xkmf(ctrl)*1.0e3,ykmf(ctrl)*1.0e3,Tb)
            plt.contour(xkmc(ctrl),ykmc(ctrl),hab(ctrl),[0])
            plt.title(r'model $|\tau _b|$ (kPa, {})'.format(year))
            plt.xticks([])
            plt.yticks([])
            plt.colorbar(pcTb)
            
            
            D_prev = D
            Tb_prev = Tb
        u_prev = u
        #Tb_prev = Tb
        
        
        
#%%
   
#files = ['ctrl.pig_sentinel_hogg.2lev.02lev.000000000000.2d.hdf5',
#         'ctrl.pig_sentinel_hogg.2lev.02lev.000000000013.2d.hdf5']   
#    
#years = [1,13]

#plot_speeds(files,years)
    
#%%    ctrl.pig_sentinel_hogg.03lev.000060000012.2d.hdf5


years = [2015,2016,2017,2018,2019,2020]
s = '{}_data/results/ctrl.pig_sentinel_{}.03lev.0000{}000016.2d.hdf5'
def fname(year,author,subdir, expt):
    s = '{}_data/{}/ctrl.pig_sentinel_{}.03lev.0000{}000016.2d.hdf5'
    return s.format(author, subdir, expt, (year-2014)*12)
    
def plot_speeds_author(author,subdir,expt):
   files=[ fname(y,author,subdir,expt) for y in years]
   plot_speeds(files,years)
   plt.savefig('pig_2016-2020-inverse-problem-{}.png'.format(expt))


#plot_speeds_author('luckman','results')
#plot_speeds_author('hogg','results')
plot_speeds_author('luckman','results_with_calving','luckman_calve')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 17:16:30 2017

@author: stephen
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
sys.path.append(os.getcwd() + '/../../../python')
from bisiclesIO import BisiclesData, write_raster, read_raster
from amrfile import io as amrio
amrio.freeAll()
slist = [(1,1,1,0.1), (1,0.5,0.5,1), (1,1,0.0,1),(0.5,1,0.5,1), (0,1,1,1) , (0,0,1,1) , (1,0,1,1), (1,0,0,1) ]
steph_lin_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('toss',slist)

def fname(year, path,  expt):
    s = '{}/plot.{}.{:06d}.2d.hdf5'
    f = s.format(path, expt, (year-2014))
    print (f)
    return f

def km(m):
    #m to km converstion
    return m*1.0e-3

def plot_hindcast_inverse(year, inverse_dir, inverse_expt, hindcast_dir, hindcast_expt):
    
        
    def read_bike(file):
        return  BisiclesData(file, level=3, croplo = (0.15,0.5), crophi = (0.35,0.75))
    
    inv = read_bike(fname(year, inverse_dir, inverse_expt))
    hin = read_bike(fname(year, hindcast_dir, hindcast_expt))# I lagged the years...
    
    fig = plt.figure(figsize=(9,3))
    
    ax = fig.add_subplot(1,3,1, aspect='equal')
    pcu = ax.pcolormesh(km(inv.x),km(inv.y) ,inv.speed, cmap=steph_lin_cmap ,vmin=0,vmax=5000)
    ax.contour(km(inv.x),km(inv.y),inv.hab,[0])
    ax.set_title('Inverse ({})'.format(year))
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(pcu)
    
    ax = fig.add_subplot(1,3,2, aspect='equal')
    pcu = ax.pcolormesh(km(inv.x),km(inv.y) ,hin.speed, cmap=steph_lin_cmap ,vmin=0,vmax=5000)
    ax.contour(km(inv.x),km(inv.y),inv.hab,[0])
    ax.set_title('Hindcast ({})'.format(year))
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(pcu)
    
    ax = fig.add_subplot(1,3,3, aspect='equal')
    pcu = ax.pcolormesh(km(inv.x),km(inv.y) ,hin.speed-inv.speed, cmap='RdBu_r' ,vmin=-500,vmax=500)
    ax.contour(km(inv.x),km(inv.y),inv.hab,[0])
    ax.set_title('Difference ({})'.format(year))
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(pcu)
    
    plt.savefig('pig_2015-2020-inverse-vs-hindcast-{}.png'.format(year),dpi=300)

years = [2015,2016,2017,2018,2019,2020]

for year in years:
    plot_hindcast_inverse(year,'luckman_data/results_with_calving','pig_sentinel_luckman_calve',
                          'hindcast_2015_2020','pig_hindcast_2015_2020')

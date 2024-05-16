#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 10:49:03 2020

@author: scornford
"""
from amrfile import io as amrio
amrio.freeAll()
import numpy as np
from netCDF4 import Dataset

import os

def _system(cmd):
    print(cmd)
    os.system(cmd)
   

def nctoamr(nc_file, hdf5_file, var_string):
    _system('nctoamr {} {} {}'.format(nc_file, hdf5_file, var_string))
    
class BisiclesData:
    """
    Store data from a BISICLES output and provide
    derived data (e.g speed from velocity)
    """
    
    def __init__(self,file_name, level=0, origin=(0,0), iord=1,
                 croplo = (0,0), crophi = (1,1), plot_file = True):
        """

        """
    
        amrid = amrio.load(file_name)
        self.time = amrio.queryTime(amrid)
        lo,hi = amrio.queryDomainCorners(amrid, level)

        lo_0 = lo
    
        for dir in [0,1]:
            L = hi[dir] - lo[dir]
            hi[dir] = int( lo_0[dir] + crophi[dir]*L )
            lo[dir] = int( lo_0[dir] + croplo[dir]*L )
            
        def read(name):
            return amrio.readBox2D(amrid, level, lo, hi, name, iord)
        
        self.x,self.y,self.thk = read('thickness')

        #adjust x,y origin
        self.x += origin[0]
        self.y += origin[1]
        
        x,y,self.usrf = read('Z_surface')
        x,y,self.topg = read('Z_base')
        
        if plot_file:
            #plot.* file data
            x,y,self.xvel = read('xVel')
            x,y,self.yvel = read('yVel')
            try:
                x,y,self.beta = read('dragCoef')
                x,y,self.hmu = read('viscosityCoef')
                x,y,self.hD = read('VIDamage')
            except:
                print ('warning: dragCoef absent, reading basal_friction instead')
                x,y,self.beta = read('basal_friction')
            
            x,y,self.acab = read('surfaceThicknessSource')
            x,y,self.ocean = read('activeBasalThicknessSource')
        else:
            #ctrl.* file data
            x,y,self.xvel = read('xVelb')
            x,y,self.yvel = read('yVelb')
            x,y,self.xvel_o = read('xVelo')
            x,y,self.yvel_o = read('yVelo')
            x,y,self.velc = read('velc')
            x,y,self.mucoef = read('muCoef')
            x,y,self.beta = read('Cwshelf')

        amrio.free(amrid)
           
        hf = np.where(self.topg < 0.0, -self.topg*1027.0/917.0, 0.0)
        self.hab = self.thk - hf
        self.speed = np.sqrt(self.xvel**2 + self.yvel**2)
        self.Tb = self.beta * self.speed
        
        if not (plot_file):
            tol = 0.1
            #self.speed_o = np.ma.masked_array( np.sqrt(self.xvel_o**2 + self.yvel_o**2), self.velc < tol)
            #gelf.misfit = np.ma.masked_array( (self.speed - self.speed_o)*self.velc, self.velc < tol)
            self.speed_o =  np.sqrt(self.xvel_o**2 + self.yvel_o**2)
            self.misfit = (self.speed - self.speed_o)*self.velc



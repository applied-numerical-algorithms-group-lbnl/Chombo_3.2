# -*- coding: utf-8 -*-
"""
A few stats to examine the output from the CalvingRate example
"""

import numpy as np
from amrfile import io as amrio

#import os.path
#os.path.exists(file_path)
class state:
    def __init__(self,x,y,usrf,bed,thk,frac,xvel,yvel,smb):
        self.x = x
        self.y = y
        self.usrf = usrf
        self.bed = bed
        self.thk = thk
        self.frac = frac
        self.xvel = xvel
        self.yvel = yvel
        self.smb = smb
        
    def u(self):
         return (self.xvel**2 + self.yvel**2)**0.5

def readplot(file,level=0):
    #read the plot file format

    s = None
    try:
        amrid = amrio.load(file)
        lo,hi = amrio.queryDomainCorners(amrid, level)
        iord = 1
        x,y,thk = amrio.readBox2D(amrid, level, lo, hi, "thickness", iord)
        x,y,frac = amrio.readBox2D(amrid, level, lo, hi, "iceFrac", iord)
        x,y,xvel = amrio.readBox2D(amrid, level, lo, hi, "xVel", iord)
        x,y,yvel = amrio.readBox2D(amrid, level, lo, hi, "yVel", iord)
        x,y,usrf = amrio.readBox2D(amrid, level, lo, hi, "Z_surface", iord)
        x,y,bed = amrio.readBox2D(amrid,level,lo,hi, 'Z_base', iord)
        x,y,smb = amrio.readBox2D(amrid, level, lo, hi, "surfaceThicknessSource", iord)
        
        amrio.free(amrid)
        s = state(x,y,usrf,bed,thk,frac,xvel,yvel,smb)
    except:
        #print ('no file ' + file)
        s = None
    return s


file_fmt = 'plot.cr_{:s}.{:s}-slip.ssa.prlim.crse1000m.{:d}lev.sg0.{:06d}.2d.hdf5'

for slip in ['no','free']:
    for lev in [0,1,2]:
        for step in range(0,1000,50):
            sx = readplot( file_fmt.format('x',slip,lev,step),lev )
            sy = readplot( file_fmt.format('y',slip,lev,step),lev )

            if (sx != None) and (sy != None) :
                delta_h = sx.thk - sy.thk.transpose()
                minv = np.min(delta_h)
                maxv = np.max(delta_h)
                l1v = np.mean(np.abs(delta_h))
                print ('amr levels: {}, step: {}, x-y h diff min = {} ,max = {}, l1 = {}'.format(lev,step,minv,maxv,l1v))
        
    

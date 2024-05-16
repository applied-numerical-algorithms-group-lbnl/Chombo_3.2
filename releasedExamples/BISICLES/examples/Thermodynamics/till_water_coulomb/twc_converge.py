from amrfile import io as amrio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import glob

amrio.freeAll()

def readplot(file, nlev):
    amrID = amrio.load(file)
    level = nlev
    #first, work out the domain corners
    lo,hi = amrio.queryDomainCorners(amrID, level)
    order = 0 # interpolation order, 0 for piecewise constant, 1 for linear
    x,y,thk = amrio.readBox2D(amrID, level, lo, hi, 'thickness', order)
    x,y,usrf = amrio.readBox2D(amrID, level, lo, hi, 'Z_surface', order)
    x,y,ux = amrio.readBox2D(amrID, level, lo, hi, 'xVel', order)
    x,y,uy = amrio.readBox2D(amrID, level, lo, hi, 'yVel', order)
    u = np.sqrt(ux**2 + uy**2) + 1.0e-10
    
    amrio.free(amrID)
    
    return x*1e-3,y*1e-3,thk,usrf,np.ma.masked_array(u,thk < 1)
    


def ctof(x):
    N = len(x)+1
    xf = np.zeros(N)
    xf[1:N] = x+0.5*(x[1] - x[0])
    return xf      

fig = plt.figure(figsize=(8,6))
n = 3
m = 2

for lev,let in zip([0,1,2,3,4],['a','b','c','d','e']):

    plt.subplot(m,n,lev+1,aspect='equal')
    file = sorted(glob.glob('plot.twc.{}lev.*.2d.hdf5'.format(lev)))[-1]
    print (file)

    x,y,h,s,u = readplot(file,lev) 

    x,y = ctof(x),ctof(y)
    #cm = plt.pcolormesh(x,y,np.log10(u),vmin = 0, vmax = 3, cmap = 'RdYlBu_r')
    cm = plt.pcolormesh(x,y,u,vmin = 1, vmax = 1000, cmap = 'RdYlBu_r', norm = LogNorm())
    plt.title(r'({}) $\Delta x =  $ {} km'.format(let,x[1] - x[0])) 
    plt.xticks([])
    plt.xlim([0,470])
    plt.ylim([0,470])
    if (lev%n == 0):
        plt.ylabel('y (km)')
        #plt.yticks([0,160,320,470])
    else:        
        plt.yticks([])

plt.subplots_adjust(wspace=0.05,hspace=0.15,right=0.95,bottom=0.05,left=0.1,top=0.95)


cbar_ax = fig.add_axes([0.7, 0.25 , 0.2 , 0.045])
cb=plt.colorbar(cm, cax=cbar_ax,   ticks=[0,1,2,3], 
                    orientation='horizontal',extend='max')  
cb.minorticks_on()
cb.set_ticklabels( [r'$10^0$',r'$10^1$',r'$10^2$',r'$10^3$'] )
cb.set_ticks([1,10,100,1000]) 
cb.set_label('Ice speed, (m/yr)')
plt.savefig('twc_converge.png',dpi=300)



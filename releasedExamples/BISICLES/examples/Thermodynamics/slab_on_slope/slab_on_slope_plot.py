from amrfile import io as amrio
import numpy as np
import matplotlib.pyplot as plt
import os

amrio.freeAll()

def readplot(file, nlayer):
    amrID = amrio.load(file)
    level = 0
    #first, work out the domain corners
    lo,hi = amrio.queryDomainCorners(amrID, level)
    hi[1] = 1
    order = 0 # interpolation order, 0 for piecewise constant, 1 for linear
    x,y,thk = amrio.readBox2D(amrID, level, lo, hi, 'thickness', order)

    sigma = np.zeros(nlayer + 2) # center of each layer, surface, base
    dsigma = 1.0 / float(nlayer) # even layers
    sigma[1:nlayer + 1] = np.arange(0.5*dsigma, 1.0, dsigma)
    sigma[nlayer + 1] = 1.0
    T = np.zeros( (len(y),len(x), len(sigma)))
    namebase = 'internalEnergy'
    for l,s in enumerate(sigma):
        
        name = namebase + '{:04d}'.format(l-1)
        if (l == 0):
            name = namebase + 'Surface'
        elif (l == nlayer + 1):
            name = namebase + 'Base'
            
        x,y,E = amrio.readBox2D(amrID, level, lo, hi, name, order)
        E = E / 2009.0
        Tpmp = 273.15 - 9.7456e-8 * 918.0 * 9.81 * thk * s
        T[:,:,l] = np.where(E > Tpmp, Tpmp, E)
    
    amrio.free(amrID)
    
    return x,y,sigma,thk,T
    
      

plt.figure(figsize=(16,8))

color = ['red','brown','purple','orange','blue','magenta','cyan','black']
#nlayer = np.array([4,8,16,32,64,128,256,512])#)128,256,512])
nlayer = np.array([8,16,32,64])
Tmin = np.zeros(len(nlayer))
nl = 0 
beta = ['0','1e6']
linestyle=['-','-.']
dx = [500,1000]

for sp,dxi in enumerate(dx):
    nti = int( 15000 / (dxi / 1000) )
    for betai,ls in zip(beta,linestyle): 
        for n,col in zip(nlayer,color):

            file = 'beta{}_{}m/plot.slab_on_slope_{}.0{}.2d.hdf5'.format(betai,dxi,n,nti)
            print (file,os.path.isfile(file))
            if (os.path.isfile(file)):
                x,y,s,h,T = readplot(file,n) 
                plt.subplot(2,1,sp+1)
                
                #TT = T[0,1,:]-273.15
                #plt.plot(TT,1-s,linestyle=ls,color='k',label=file)
                
                TT = T[0,0,:]-273.15
                plt.plot(TT,1-s,linestyle=ls,color=col,label=file)
               

def annotate(sp):
    plt.subplot(2,1,sp)  
    plt.legend()   

    plt.xlabel(r'$T$($^\circ$C)')
    plt.ylabel(r'$1- \sigma$')
    plt.ylim(0,0.4)
    plt.xlim(-37,-27)
    
    
annotate(1)
annotate(2)

plt.savefig('T_sigma.png',dpi=300)



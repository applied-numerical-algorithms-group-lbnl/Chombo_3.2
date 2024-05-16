from amrfile import io as amrio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.patches as pat


plt.figure(figsize=(4,6))
plt.subplot(111,aspect='equal')

hmin = 0.0
hmax = 3000.0

amrio.freeAll()
amrID = amrio.load("plot.amundsen.2d.hdf5")
time = amrio.queryTime(amrID)



n_lev = amrio.queryLevelNumber(amrID)
lev_col = ['grey','blue','purple','red','orange']

for lev in range(0,n_lev):
    n_fab = amrio.queryFABNumber(amrID,lev)
    print(lev, n_fab)
    for fab in range(0,n_fab):
        x,y,h = amrio.readFAB(amrID,lev,fab,0)
        dx = (x[1]-x[0])
    
    
    
        x0 = x[0] - 0.5*dx
        y0 = y[0] - 0.5*dx
        x1 = x[-1] + 0.5*dx
        y1 = y[-1] + 0.5*dx
    
        eps = dx * 0.1
        xx = np.arange(x0,x1+eps,dx)
        yy = np.arange(y0,y1+eps,dx)
        plt.pcolormesh(xx,yy,h,vmin=hmin,vmax=hmax)
    
        plt.plot( [x0,x0,x1,x1,x0],[y0,y1,y1,y0,y0], color=lev_col[lev], 
                 lw=0.5, label = r'$\Delta x = ${} m'.format(dx))

plt.legend()

amrio.free(amrID)
plt.savefig("libamrfile_python_mesh.png")

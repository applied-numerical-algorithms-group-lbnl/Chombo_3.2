from amrfile import io as amrio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.patches as pat

amrID = amrio.load("plot.amundsen.2d.hdf5")

thkcomp = "thickness"
thklim = col.Normalize(0.0,4000.0) # limits for thickness colormap
thkc = [0,1000,1500,2000]
#read a box of thickness data at the lowest resolution
level = 0
#first, work out the domain corners
lo,hi = amrio.queryDomainCorners(amrID, level)
order = 0 # interpolation order, 0 for piecewise constant, 1 for linear
x0,y0,thk0 = amrio.readBox2D(amrID, level, lo, hi, thkcomp, order)

#set up figure axes
asp = (max(y0)-min(y0))/(max(x0)-min(x0))
fig = plt.figure(1,figsize=(6, 6*asp))
plt.xlim (min(x0),max(x0))
plt.ylim (min(y0),max(y0))
plt.xticks([0,250e+3,500e+3])
plt.yticks([0,250e+3,500e+3],rotation=90)

#color and contour plot
plt.pcolormesh(x0,y0,thk0,norm=thklim,figure=fig)
cs = plt.contour(x0,y0,thk0,thkc,figure=fig,norm=thklim,colors='black')

plt.clabel(cs, inline=1, fontsize=10)
#read thickness data at level 1 resolution
lo = [50,50]
hi = [150,150]
level = 1
x1,y1,thk1 = amrio.readBox2D(amrID, level, lo, hi, thkcomp, order)
plt.pcolormesh(x1,y1,thk1,figure=fig,norm=thklim)
plt.contour(x1,y1,thk1,thkc,figure=fig,norm=thklim,colors='black')

#rectangle around the highres area
dx = x1[1] - x1[0]
c=[min(x1)-dx/2.0,min(y1)-dx/2.0]
w = max(x1)-min(x1) + dx/2.0
h = max(y1)-min(y1) + dx/2.0
plt.gca().add_patch(pat.Rectangle((min(x1)-dx/2.0,min(y1)-dx/2.0) , w, h, edgecolor = 'pink', fill=False))




plt.savefig("libamrfile_python.png")

time = amrio.queryTime(amrID)
#print "time = ", time

amrio.free(amrID)

import math


ylen = 80.0e+3
xlen = 80.0e+3

def thickness(x,y):

    x0 = x
    y0 = y
    w = width(x)
    x = x/1.0e+3
    y = (y-ylen/2.0)/1.0e+3


    thickness = 0.0
    surface = 0.0
    
    smax = 1.0*(84.3 * (math.sqrt (69-0.0)))

    if ((x <= 69.0)  & (topography(x0,y0) < 1000.0)):
        surface = min(smax,1.0*(84.3 * (math.sqrt (69-x))))

    if surface > 0.0:
        thickness = surface - topography(x0,y0)
        

    return thickness 

def width(x):

    x = x/1.0e+3

    xx = max(0.0,x)
    w = max(4.0, ((2*1.236e+4) / math.pow((8.964 + xx),2.629)))
 
    return w 


def inside(x,y):
# return True if (x,y) is inisde the channel, False otherwise
    r = False
    hw = 1.0e+3 * width(x) * 0.5
    if (abs(y-0.5*ylen) < hw):
        r = True

    return r

def outside(x,y):
    return ( not (inside(x,y)) )

def tag(x,y,dx,thck,topg):

    tag = -1.0
    y= y - ylen/2.0

    #limit finest refinement to the downstream region
    if (x > 24.0e+3):
    #refine the whole of the 4 km stripe to 125m
        w= 4.0e+3       
        if ( (dx > 250.0) & ((abs(y) - dx/2.5) < w/2.0) & (thck > 0.0)):
            tag = 1.0     
    #refine the shear margins up to the maximum
        if ( (  thck > 0.0 ) & ( w/2.0 - abs(y) < 2.0*dx ) ):
            tag = 1.0
    #compute the thickness above flotation and refine the ice shelf
    #and regions close to the gl up to the maximum
        hab = min(thck,thck+topg*1028/918)
        if ( (thck > 1.0) & (hab < 100.0)):
            tag = 1.0


    #refine the curve to 500m
    hw = 1.0e+3*width(x)/2.0
    hy = abs(y)
    if ( (dx > 500.0) & ( abs(hy - hw) < (dx*2.5) ) ):
        tag = 1.0

    #refine the majority of the ice filled domain to 500m    
    #if ( (dx > 500.0) & ((x > 16.0e+3) & (  thck > 0.0 ))):
    #    tag = 1.0

    return tag

def topography(x,y):
    
 
    pi = 3.141592653
    w = width(x)
    x = x/1.0e+3 
    y = (y-ylen/2.0)/1.0e+3

    topography = 0.0

    if (abs(y) < w/2.0):
        topography = ((-550.0/80.0)*x + 50.0)
    else:
        topography = 2000.0

    return topography
            


def accumulation(x,y,t,thck,topg):
    acc = 0.0
    x = x/1.0e+3 
    y = (y-ylen/2.0)/1.0e+3
     
    acc = 0.0

    if topg < 1000.0:
        acc = max(-3.5, (x*(-8.0/43.0) + 3.0))
        
    if (x > 71.0e+3):
        acc = -10000.0
        

    return acc




def rhs(x,y,dx,idir,rhs,rhogH,s,sleft,sright):
#modify the rhs = rho*g*H*ds/dx[dir] to avoid
#large slopes at stream / wall interfaces. The
#default implementation breaks down when a coarse-fine
#interface intersects the wall.     
   # print "testing rhs(",dx,x,y,") "


    xleft = x
    xright = x
    if (int(idir) == 0):
        xleft = x - dx
        xright = x + dx
    
    yleft = y
    yright = y
    if (int(idir) == 1):
        yleft = y - dx
        yright = y + dx


    if (inside(x,y)):
        alter = False

        if ( not (inside(xleft,yleft))):
            sleft = min(sleft,s)
            alter = True

        if ( not (inside(xright,yright))):
            sright = min(sright,s)
            alter = True
            
        if (alter):
            #print "alter rhs(",dx,x,y,") "
            rhs = 0.5 * rhogH * (sright-sleft) / dx
    else:
        rhs = 0.0

    return rhs

def facevel(x,y,dx,idir,facevel):
    #set facevel to zero outside the ice stream and at its edges
    #x and y will be face centered
    hdx =  0.5 * dx

    if (int(idir) == 0):
        if ( outside(x-hdx,y) or outside(x+hdx,y) ):
            facevel = 0.0

    if (int(idir) == 1):
         if ( outside(x,y-hdx) or outside(x,y+hdx) ):
             facevel = 0.0       
   

    return facevel

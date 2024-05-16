import math


ylen = 40.0e+3
xlen = 80.0e+3

def thickness(x,y):

    x0 = x
    y0 = y
    w = width(x)
    x = x/1.0e+3
    y = (y-ylen/2.0)/1.0e+3


    thickness = 0.0
    surface = 0.0

    if (topography(x0,y0) < 1000.0):
    
        if (x <= 45.0):
            surface = 1.0*(84.3 * (math.sqrt (45-x)))
            thickness = surface - topography(x0,y0)

        elif (x <= 55.0):
            thickness = 100.0

        elif (x <= 65.0):
            surface = 0.0   

        elif (x <= 75.0):  
            thickness = 100.0  

    return thickness 

def width(x):

    x = x/1.0e+3
    w = max(2.0, (1.236e+4 / math.pow((8.964 + x),2.629)))

    return w 

def tag(x,y,dx,thck,topg):

    tag = -1.0

    #refine the whole of the 2 km stripe to 250m
    w= 2.0e+3
    y= y - ylen/2.0

    if ( (dx > 250.0) & ((abs(y) - dx/2.0) < w/2.0) & (thck > 0.0)):
        tag = 1.0     

    #refine the curve to 500m
    hw = 1.0e+3*width(x)/2.0
    hy = abs(y)
    if ( (dx > 500.0) & ( abs(hy - hw) < dx ) & (  thck > 0.0)):
        tag = 1.0


    #compute the thickness above flotation and refine the ice shelf
    #and regions close to the gl up to the maximum
    hab = min(thck,thck+topg*1028/918)
    if ( (thck > 1.0) & (hab < 100.0)):
        tag = 1.0


    return tag

def topography(x,y):
    
 
    pi = 3.141592653
    w = width(x)
    x = x/1.0e+3 
    y = (y-ylen/2.0)/1.0e+3

    topography = 0.0

    if (abs(y) < w/2.0):
        if (x < 20.0):
            topography = ((-15.0/2.0) * x) + 50.0
        elif (x < 44.0): 
            topography = -250.0 + (150.0*math.cos(2.0*pi*(((x-20)/24.0))))
        elif (x  < 68.0): 
            topography = -250.0 + (150.0*math.cos(2.0*pi*(((x-44)/24.0))))
        else: 
            topography = -250.0 + (150.0*math.cos(2.0*pi*(((x-68)/24.0))))
    else:
        topography = 9500.0

    return topography
            


def accumulation(x,y,t,thck,topg):
    acc = 0.0
    x = x/1.0e+3 
    y = (y-ylen/2.0)/1.0e+3
     
    acc = 0.0

    if topg < 1000.0:
        acc = max(-4.2, (x*(-8.0/44.0) + 3.0))
        
    if (x > 71.0e+3):
        acc = -10000.0
        

    return acc

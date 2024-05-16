import math

def slope(x):
    return 50.0 - 4.0 * x / 1.0e+3

def width(x):
    return 4.0e+3 + 16.0e+3 * math.exp(-x/16.0e+3)

def tag(x,y,dx,thck,topg):
    w = width(x) + 1.0e+2
    y = y - 32.0e+3
    tag = -1.0 #don't tag unless...
    if ((dx > 250.0) & ((abs(y) - dx/2) < w)):
        tag = 1.0
        
    return tag

def thickness(x,y):
    thickness = 0.0
  
    if (x < 60.0e+3):
        surface = 300.0 + slope(x)
        thickness = surface - topography(x,y)
        thickness = max(0.0, thickness)


        thickness = 200.0
    return thickness

def topography(x,y):
    w = width(x)
    y = y - 32.0e+3
    topography = slope(x) # + 1.0e+3 - 1.0e+3 * math.exp(-(y/w * y/w))
    
    return topography


def accumulation(x,y,t,thck,topg):
    acc = -1000.0
    surface = 300.0 + slope(x)
    if ((x < 80.0e+3) & (surface > topg)):
        acc = 0.5 - 0.5 * x / 70e+3
    
    return acc

def friction(x,y,t,thck,topg):
    friction = 1.0e+6
    surface = 300.0 + slope(x)
    if (  surface > topg ):
        friction = 1.0e+3
        
    return friction

def frictionTemp(x,y,t,thck,topg):
    friction = 1.0e+6
    surface = 300.0 + slope(x)
    if (  surface > topg ):
        friction = 1.0e+4
        
    return friction

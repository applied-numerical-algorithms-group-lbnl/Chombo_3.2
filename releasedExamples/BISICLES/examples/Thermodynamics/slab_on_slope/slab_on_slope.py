#constant thickness slab on slope
L = 64.0e+3
ACAB = 4.0
THICKNESS = 1000.0
ROS = ACAB / THICKNESS * 1.1


import random
random.seed(0)

def thickness(x,y):
    return 1000.0

def topography(x,y):
    return 1000.0 * (1.0 - x/L)

def clean_velocity(x,y):
    return ROS * x,0.0

def noisy_velocity(x,y):
    r = NOISE * (2.0*random.random() - 1.0) * x
    u,v = clean_velocity(x,y)
    return u + r , v


def velocity(x,y,*etc):
    return clean_velocity(x,y)

def stemp(x,y,t,thck,topg,*etc):
    s = thck + topg
    return 273.15 - 35.0

def bflux(x,y,t,thck,topg,*etc):

    #return 100e-3 * 3600 * 24 * 365.24
    return 150e-3 * 3600 * 24 * 365.24

def acab(x,y,t,*etc):
    return ACAB

import math

def frictioncoef(x,y,*etc):

    friction = 1.0e+6
    if (abs(y-64.0e+3) < 8.0e+3):
        friction = 0.0e+0
    return friction

def frictioncoef_spot(x,y,*etc):
    friction = 1.0e+6
    if (abs(y-64.0e+3) < 8.0e+3):
        friction = 0.0e+0

    if (math.sqrt(math.pow(y-64.0e+3,2) + math.pow(x-64.0e+3,2)) < 4.0e+3):
        friction = 100.0

    return friction


def temperature(x,y,thck,topg,sigma):
    T = 263.0
    if (abs(y-64.0e+3) < 8.0e+3):
        T = T + sigma*10.0
    return T

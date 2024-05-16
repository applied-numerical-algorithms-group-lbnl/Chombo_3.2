# LxL periodic domain flat bed, sinusoidal thickness
L = 64.0e+3
import math
k = 2 * math.pi / L

def thickness(x,y):
    return 1000.0 + 500.0 * math.sin(k * x)

def topography(x,y):
    return 1000.0

def velocity(*etc):
    return 100.0

def velocityB(x,y,*etc):
    return 100.0 + 50.0 * math.sin(k * x)

def thicknessB(x,y,*etc):
    return 1000.0 


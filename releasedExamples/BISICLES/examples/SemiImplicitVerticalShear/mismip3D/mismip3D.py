#MISMIP 3D geoemtry imposed through the python interface
import math

def slope(x):
    return -100.0 - x * 1.0e-3

#spin up from 150 m
def thickness(x,y):
    return 150.0

def topography(x,y): 
    return slope(x)


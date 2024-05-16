import math as m

Ldomain = 480e+3
dxCrse = 10.0e+3

def basal_flux(x,y,*etc):
    #maintain calving front at the bottom and right edges,
    #aligned with a cell face on the coarsesst grid
    huge = 1.0e+4
    f = 0.0
    
    if (x > Ldomain - dxCrse):
        f = -huge

    if (y > Ldomain - dxCrse):
        f = -huge
        
    return f
def basal_flux_reverse(x,y,*etc):
    return -basal_flux(x,y,*etc)

def thickness(x,y,*etc):
    #initial thickness
    return 100.0

def topography(x,y,*etc):
    return 1.0


def till_water_drain_factor(x,y,*etc):

    factor = 0.005
    fast_factor = 1.0e+3
    r2 = x**2 + y**2
    if (r2 < 0.25 * Ldomain**2):
        factor = fast_factor

    return factor

import math as m

#shallow 
B0 = -150.0 #vs 0 from gael
B2 = -728.8 #from gael
B4 = 343.91 #from gael
B6 = - 50.57 #from gael


def Bx(x):
    xx = x/300.0e3
    xx2 = xx*xx
    xx4 = xx2*xx2
    xx6 = xx4*xx2
    return B0 + B2*xx2 + B4*xx4 + B6*xx6 

def By(y,fc,dc,wc):
    return dc/(1.0+m.exp(-2.0*(y-wc)/fc)) + dc/(1.0+m.exp(2.0*(y+wc)/fc))

def slope(x):
    return 50.0 - 4.0 * x / 1.0e+3

def width(x):
    return 4.0e+3 + 16.0e+3 * math.exp(-x/16.0e+3)

def tag(x,y,dx,thck,topg):
    tag = -1.0 #don't tag unless...
    return tag

def thickness24(x,y):
    thickness = 0.0
  
    if (x < 640.0e+3):
        thickness =  100.0
        
    return thickness


def thicknessinf(x,y):
    thickness = 0.0
  
    if (x < 640.0e+3):
        thickness =  1000.0
        
    return thickness



#infinite with topography, ie 1D problem
def topographyinf(x,y):
    topography = Bx(x) 
    topography = max(topography, -720.0)
    return topography

#24 km channel half-width 
def topography24(x,y):
    y = y - 40e+3
    topography = Bx(x) + By(y,4.0e+3,5.0e+2,24.0e+3)
    topography = max(topography, -720.0)
    return topography

#16 km channel half-width 
def topography16(x,y):
    y = y - 40e+3
    topography = Bx(x) + By(y,4.0e+3,5.0e+2,16.0e+3)
    topography = max(topography, -720.0)                  
    return topography 

#12 km channel half-width 
def topography12(x,y):
    y = y - 40e+3
    topography = Bx(x) + By(y,4.0e+3,5.0e+2,12.0e+3)
    topography = max(topography, -720.0) 
    return topography 


def constfriction_24(x,y,t,thck,topg):
    friction = 1.0e4
    return friction

def constfriction_16(x,y,t,thck,topg):
    friction = 1.0e4
    return friction


def constfriction_12(x,y,t,thck,topg):
    friction = 1.0e4
    return friction

def constfriction_inf(x,y,t,thck,topg):
    friction = 1.0e4
    return friction

def hardfriction_24(x,y,t,thck,topg):
    friction = 1.0e5
    return friction

def hardfriction_16(x,y,t,thck,topg):
    friction = 1.0e5
    return friction


def hardfriction_12(x,y,t,thck,topg):
    friction = 1.0e5
    return friction

def hardfriction_inf(x,y,t,thck,topg):
    friction = 1.0e5
    return friction




def chanfriction_24(x,y,t,thck,topg):
    y = y - 40e+3
    friction = (2.0e4) *  By(y,4.0e+3,5.0e+2,12.0e+3)/By(40e+3,4.0e+3,5.0e+2,24.0e+3) + 0.5e4
    return friction

def chanfriction2_24(x,y,t,thck,topg):
    y = y - 40e+3
    friction = (1.75e4) *  By(y,4.0e+3,5.0e+2,12.0e+3)/By(40e+3,4.0e+3,5.0e+2,24.0e+3) + 0.75e4
    return friction

def chanfriction_16(x,y,t,thck,topg):
    y = y - 40e+3
    friction = (2.0e4) *  By(y,4.0e+3,5.0e+2,12.0e+3)/By(40e+3,4.0e+3,5.0e+2,16.0e+3) + 0.5e4
    return friction

def chanfriction2_16(x,y,t,thck,topg):
    y = y - 40e+3
    friction = (1.75e4) *  By(y,4.0e+3,5.0e+2,12.0e+3)/By(40e+3,4.0e+3,5.0e+2,16.0e+3) + 0.75e4
    return friction



def chanfriction_12(x,y,t,thck,topg):
    y = y - 40e+3
    friction = (2.0e4) *  By(y,4.0e+3,5.0e+2,12.0e+3)/By(40e+3,4.0e+3,5.0e+2,12.0e+3) + 0.5e4
    return friction

def chanfriction2_12(x,y,t,thck,topg):
    y = y - 40e+3
    friction = (1.75e4) *  By(y,4.0e+3,5.0e+2,12.0e+3)/By(40e+3,4.0e+3,5.0e+2,12.0e+3) + 0.75e4
    return friction

#big melr rate to hold calving front at x = 640 km
def cfsource(x):
    m = 0.0
    if (x > 640.0e+3):
        m = 1.0e+4
    return -m

def spinbasalsource(x,y,*etc):
    melt = 0.0
    return -melt + cfsource(x)



def melt1basalsource(x,y,t,thck,topg):
#m = -5.0e-2 * z_bottom * (T-Tf) * tanh( (z_bottom - z_bathymetry)/ 75)                                                                      
    rhoi = 918.0
    rhoo = 1028.0
    Ta = 1.0
    Cz =  7.61e-4 
    mr = 0.0
    zb = -thck * rhoi / rhoo
    COEF=5.0e-2

    if ( (zb < 0.0) & (t < 210.0) & (t >= 10.0)):
        wct = zb-topg                                                                       
        Tf =  -1.85 + Cz*zb
        mr = COEF * zb * (Tf - Ta) * m.tanh(wct/75.0) 

    return -mr + cfsource(x)

def melt2basalsource(x,y,t,thck,topg):
#m = -5.0e-2 * (z_bottom+100) * (T-Tf) * tanh( (z_bottom - z_bathymetry)/ 75)                                                                      
#not melting above -100 m
    rhoi = 918.0
    rhoo = 1028.0
    Ta = 1.0
    Cz =  7.61e-4 
    mr = 0.0
    zb = -thck * rhoi / rhoo
    COEF=5.0e-2

    if ( (zb < -100.0) & (t < 210.0) & (t >= 10.0)):
        wct = zb-topg                                                                       
        Tf =  -1.85 + Cz*zb
        mr = COEF *  (zb+100.0) * (Tf - Ta) * m.tanh(wct/75.0) 

    return -mr + cfsource(x)

def melt3basalsource(x,y,t,thck,topg):
# simplified melt, no (weak) T-Tf dependency
#m = Omega * (z_bottom - z_0) * tanh( (z_bottom - z_bathymetry)/ 75)            
#not melting above z_0
    z_0 = -100.0
    rhoi = 918.0
    rhoo = 1028.0
    zb = -thck * rhoi / rhoo
    Omega = -0.2
    mr = 0.0
    if ( (zb < z_0) & (t < 210.0) & (t >= 10.0)):
        wct = zb-topg                                                           
        mr = - Omega * (zb - z_0) * m.tanh(wct/75.0) 

    return mr + cfsource(x)

def melt4basalsource(x,y,t,thck,topg):
# simplified melt, no (weak) T-Tf dependency
#m = Omega * (z_bottom - z_0) * tanh( (z_bottom - z_bathymetry)/ 75)            
#not melting above z_0
#same as melt 3, but turn off after 100 years
    z_0 = -100.0
    rhoi = 918.0
    rhoo = 1028.0
    zb = -thck * rhoi / rhoo
    Omega = -0.2
    mr = 0.0
    if ( (zb < z_0) & (t < 110.0) & (t >= 10.0)):
        wct = zb-topg                                                           
        mr = - Omega * (zb - z_0) * m.tanh(wct/75.0) 

    return mr + cfsource(x)


def melt5basalsource(x,y,t,thck,topg):
# melt-rate that emulates a calving event
# at x = 480 km by imposing 100 m/a melt
    m = cfsource(x)
    if ( ( t < 110.0) & (t >= 10.0) ):
        if (x > 480.0e+3):
            m = -1.0e+2

    return m 

def melt42basalsource(x,y,t,thck,topg):
# simplified melt, no (weak) T-Tf dependency
#m = Omega * (z_bottom - z_0) * tanh( (z_bottom - z_bathymetry)/ 75)            
#not melting above z_0
#same as melt 4, but don't turn off
    z_0 = -100.0
    rhoi = 918.0
    rhoo = 1028.0
    zb = -thck * rhoi / rhoo
    Omega = -0.2
    mr = 0.0
    if ( (zb < z_0) & (t >= 10.0)):
        wct = zb-topg                                                           
        mr = - Omega * (zb - z_0) * m.tanh(wct/75.0) 

    return mr + cfsource(x)


def melt52basalsource(x,y,t,thck,topg):
# melt-rate that emulates a calving event
# at x = 480 km by imposing 100 m/a melt
# same as melt 5, but don't turn off
    m = cfsource(x)
    if ( t >= 10.0 ):
        if (x > 480.0e+3):
            m = -1.0e+2

    return m 
    
def melt0basalsource(x,y,t,thck,topg):
# just the spinu source, for control runs
    return spinbasalsource(x,y,t,thck,topg)




def meltxylarbasalsource(x,y,t,thck,topg):
#melt4 plus Xylar's calving front, delayed by 100 years
    t = t - 90.0
    mr = melt4basalsource(x,y,t,thck,topg)

    #now use Xylar's front
    mcr = 0.0
    y = y - 40.0e+3
    step = 0.0
    if (abs(y) < 32.0e+3):
        step = 1.0

    xcy = 548.0e+3 + 92.0e+3 * 0.5 * step
 
    if (x > xcy):
        mcr = -1.0e+4

    return mr + mcr

def melthlimbasalsource(x,y,t,thck,topg):
#melt4 plus removal of ice < 100m thick
    t = t - 90.0
    mr = melt4basalsource(x,y,t,thck,topg)

    #remove thin ice
    mcr = 0.0
    if (thck < 100.0):
        mcr = -1.0e+4

    return mr + mcr

def melthlim80basalsource(x,y,t,thck,topg):
#melt4 plus removal of ice < 80m thick
    t = t - 90.0
    mr = melt4basalsource(x,y,t,thck,topg)

    #remove thin ice
    mcr = 0.0
    if (thck < 80.0):
        mcr = -1.0e+4

    return mr + mcr

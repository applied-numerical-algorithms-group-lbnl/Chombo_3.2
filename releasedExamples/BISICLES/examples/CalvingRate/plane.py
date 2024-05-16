from math import tanh,pi

L = 128.0e+3 #domain length
W = 8.0e+3 #domain with

b0 = -128.0 # bedrock elevation at divide
dbdx = -1.0 / 256 * 4.0  # bedrock slope in x

a0 = 4.0 # accumulation at divide
dadx = -8.0 / L # accumulation slope

cr0 = 0.0 # calving rate at divide
crL = 512.0 # calving rate at L  (no slip version)
crL = 1024.10 # calving rate at L (free slip version)
crL = 768.0

dcdx = (crL - cr0) / L

T_SPIN = -1.0
T_SWITCH = 15000.0

u0 = 0.0
dudx = 0.75*dcdx

def velocity_x_free_slip(x,y,*etc):
    xstar = 2.0 * pi * x / L 
    u = 0.5 * crL * tanh(xstar)
    return (u,0.0)

def velocity_x_no_slip(x,y,*etc): 
    u,v = velocity_x_free_slip(x,y,*etc)
    ystar = 2.0 * y/W -1.0
    #u = 0.5 * (u + u * (1.0 - ystar*ystar))
    u = u * (1.0 - ystar*ystar)
    v = - ystar * max(0.0, acab_x(x,y,*etc))
    return (u,v)

def velocity_y_free_slip(x,y,*etc):
    return (0.0,velocity_x_free_slip(y,x,*etc)[0])

def velocity_y_no_slip(x,y,*etc):
    return (0.0,velocity_x_no_slip(y,x,*etc)[0])

def thickness(x,y):
    h = 0.0
    if (x < 500.0e+3):
        h = 256.0
    
    return h

def topography_x(x,y):
    return b0 + dbdx * x

def topography_y(x,y):
    return topography_x(y,x,*etc)

def constfriction(x,y,t,thck,topg):
    return 1.0e+4

def acab_x(x,y,*etc):
    return (a0 + dadx * x)

def acab_y(x,y,*etc):
    return acab_x(y,x,*etc)

def calving_rate_x(x,y,t,*etc):
    cr = 0.0
    if (t > T_SPIN):
        cr = cr0 + dcdx *x 
    return  cr 

def calving_rate_ux(x,y,t,*etc):
    cr = 0.0
    if (t > T_SPIN):
        cr = (2.0 * x/L)**3
    if (t > T_SWITCH):
        cr = (3.0/2.0 * x/L)**3
        
    return  cr 

def calving_rate_y(x,y,t,*etc):
    return  calving_rate_x(y,x,t,*etc) 

def calving_rate_uy(x,y,t,*etc):
    return  calving_rate_ux(y,x,t,*etc)

def spinbasalsource(x,y,*etc):
    return 0.0


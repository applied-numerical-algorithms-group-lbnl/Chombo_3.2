import math as m

L = 64.0e+3
W = 128.0e+3

L1 = 0.75 * L #gl to crack
L2 = 2.0e+3 #crack width
L3 = L - L1 - L2 #crack to front

W1 = 0.25 * W
W2 = 0.5  * W
W3 = W - W1 - W2

#pinning point 
r_pin = 4.0e+3
x_pin = W - W3 + r_pin
y_pin = L1 
rsq_pin  = r_pin**2

#berg in end of rift
r_berg = r_pin
x_berg = W - W3
y_berg = L1 
rsq_berg  = r_berg**2


def topography(x,y):

    b = -1000.0

    rsq = (x - x_pin)**2 + (y - y_pin) **2 

    if rsq < rsq_pin :
        b = -200.0
    
    return b

def thickness(x,y):

    h = 600.0
    h = h + 16.0 * ( (x-W/2.0) / W )**2.0
    h = h - 300 * y/L
    
    if ( y > L):
        h = 0.0
    
    if ((x > W1) and (x < W - W3)):
        if (( y > L1) and (y < L - L3)):
            h = 0.0
    elif (x > W2):
        if (y > L1):
            h = 0.0

    rsq = (x - x_berg)**2 + (y - y_berg) **2         
    if (rsq < rsq_berg):
        h = 100.0
    
            
    return h
    

def beta(x,y,*etc):
    return 1.0e+3

def no_thin(x,y,t,thck,topg):

    m = 0.0
    if (thck < -10.0):
        m = -100.0

    return m

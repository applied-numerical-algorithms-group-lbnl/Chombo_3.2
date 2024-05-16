def sum(a,b):
    return a+b

def minus(a,b):
    return a-b

def mult(a,b):
    return a*b

def surfaceAbove(usrf,thck):
    #surface above flotation (assuming densities)
    return usrf,thck,usrf - (1.0-918.0/1028.0)*thck

def mod2D(u,v):
    return pow(u*u + v*v,0.5)

def mod3D(u,v,w):
    return pow(u*u + v*v + w*w,0.5)


def linearToPlastic(C,u,v):
    # solve C|u| = B|u|^{1/3} for B 
    umod = mod2D(u,v)
    return C  * (1+pow(umod, 2.0/3.0)), C, C*umod

def plasticToLinear(B,u,v):
    # solve C|u| = B|u|^{1/3} for C
    umod = mod2D(u,v)
    return B  * pow(umod+1, -2.0/3.0)



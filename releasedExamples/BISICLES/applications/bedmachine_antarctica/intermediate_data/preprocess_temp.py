import numpy as np
from netCDF4 import Dataset
    
def sigma(n_layer):
    n = n_layer+1 # no of faces
    tol = 1.0e-10
    x = np.arange(0,1+tol,1.0/n)
    f = np.pi/2.0
    s_face = np.tanh(f*x)/np.tanh(f)
    s_mid = 0.5*(s_face[0:n-1] + s_face[1:n])
    return s_face, s_mid
    

def morlighem_temperature_nc(out_temperature_file, geometry_file, morlighem_temperature_file, sigma_mid, T_MAX=273.0, T_INCR=0.0):

    from  scipy.interpolate import RectBivariateSpline
    from scipy import ndimage
    
    ncgeo = Dataset(geometry_file,'r')
    x = ncgeo.variables['x'][:]
    y = ncgeo.variables['y'][:] 
    thk = ncgeo.variables['thk'][:]
    topg = ncgeo.variables['topg'][:]
    ncgeo.close()

    ncmor =  Dataset(morlighem_temperature_file,'r')
    xm = ncmor.variables['x'][:]
    ym = ncmor.variables['y'][:]
    zeta = ncmor.variables['zeta'][:]
    sm = 1.0-np.flipud(zeta)
    nmor = len(sm)
    
    ncout = Dataset(out_temperature_file,'w')
    xdim = ncout.createDimension('x',size=len(x))
    ydim = ncout.createDimension('y',size=len(y))
    xv = ncout.createVariable('x','f8',('x'))
    yv = ncout.createVariable('y','f8',('y'))
    xv[:] = x
    yv[:] = y
    
    rhoo = 1027.0
    rhoi = 917.0
    s = thk + topg
    sf = (1.0 - rhoi/rhoo) * thk
    s = np.where( s > sf, s, sf) 


    ice = np.where(thk < 1, 0, 1)
    dx = x[1]-x[0]
    ice_not_edge = ndimage.minimum_filter(ice,int(8e3/dx))
    #thkm = RectBivariateSpline(x,y,thk,kx=1,ky=1)(xm,ym)

    

    def read_morlighem_layer(x,y,k):
        layer = nmor - k  -1
        print('reading layer {}, zeta = {}, sigma = {}'.format(layer, zeta[layer], 1-zeta[layer]))
        T0 = ncmor.variables['temperature'][layer,:,:]
        T0 = np.where(T0 < T_MAX, T0, T_MAX)
        #T0 = np.where(thkm > 1 , T0, T_MAX)
        print ('interpolating x,y')
        T = RectBivariateSpline(xm,ym,T0,kx=1,ky=1)
        Txy = T(x,y) + T_INCR
        Txy  = np.where(ice_not_edge, Txy, T_MAX)
        print ('    ...filtering')

        Txyf = ndimage.maximum_filter(Txy, int(12e3/dx))
        #Txyf = ndimage.gaussian_filter(Txyf, int(16e3/dx))
        Txy = np.where(ice_not_edge < 1, Txyf, Txy)
        
        return np.where(thk > 1, Txy, T_MAX)

    kp,kp_prev = 1,0
    Tm,Tp = None,None
    for layer in range(0, len(sigma_mid)):

        s = sigma_mid[layer]
        while s > sm[kp]:
            kp += 1

        ds = sm[kp]-sm[kp-1]
        wm = (sm[kp]-s)/ds
        wp = 1.0 - wm
        print ('sigma = {}, kp = {}, sm[kp-1] = {}, sm[kp] = {}, wm = {}'.format(s,kp,sm[kp-1],sm[kp],wm)) 
            
        if (kp != kp_prev):
            kp_prev = kp
            if (type(Tp) != type(None)):
                Tm = Tp.copy() # i should swap pointers. too lazy
            else:
                Tm = read_morlighem_layer(x,y,kp-1)  
            Tp = read_morlighem_layer(x,y,kp)

            print (type(Tp),type(Tm))
            
       
        
        T = wm*Tm + wp*Tp

        #T = np.where(T < T_MAX, T, T_MAX)
        
        name = 'temp{:06d}'.format(layer)
        Tv = ncout.createVariable(name,'f8',('y','x'))
        Tv[:,:] = T
    
    ncout.close()

    
    

def simple_temperature_nc(temperature_file, geometry_file, sigma_mid):
    

    T_0 = 260
    T_MAX = 273
    LAPSE_RATE = -8.0e-3
      
    ncgeo = Dataset(geometry_file,'r')
    x = ncgeo.variables['x'][:]
    y = ncgeo.variables['y'][:] 
    thk = ncgeo.variables['thk'][:]
    topg = ncgeo.variables['topg'][:]
    ncgeo.close()
    
    ncout = Dataset(temperature_file,'w')
    xdim = ncout.createDimension('x',size=len(x))
    ydim = ncout.createDimension('y',size=len(y))
    xv = ncout.createVariable('x','f8',('x'))
    yv = ncout.createVariable('y','f8',('y'))
    xv[:] = x
    yv[:] = y
    
    rhoo = 1027.0
    rhoi = 917.0
    s = thk + topg
    sf = (1.0 - rhoi/rhoo) * thk
    s = np.where( s > sf, s, sf) 
    
    for layer in range(0, len(sigma_mid)):
        
        def z(s):
            return s - 0.5* ( sigma_mid[layer] )*thk
        
        T = T_0 + LAPSE_RATE * z(s)
        T = np.where(T > T_MAX, T_MAX, T)
        
        name = 'temp{:06d}'.format(layer)
        Tv = ncout.createVariable(name,'f8',('y','x'))
        Tv[:,:] = T
    
    ncout.close()


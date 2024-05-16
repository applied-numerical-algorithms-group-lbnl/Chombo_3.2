### Inverse problem post-processing.
### The inverse problem finds c_one, which
### is just telling us the basal friction
### |Tb| = c_one | u | at one point in time
### But we don't normally want to do prognostics runs
### with that rule, we want another rule e.g 
### |Tb| = c_third | u |^1/3
### that produces the same |Tb|
### so here, we work out c_third from |u| and c_one.
###
### To use: define $extract, $flatten, and $pythonf (filetools) and the
### solution to the inverse problem $soln (an hdf5 file e.g inverse_0/regC15e3mu15e10/ctrl.ase_bmach.03lev.000000000055.2d.hdf5)
### then run
### > $extract $soln tmp.2d.hdf5 Cwshelf muCoef xVelb yVelb xVelo yVelo velc
### > $flatten tmp.2d.hdf5 flat.2d.hdf5 3 # level 3: 500m
### > $pythonf flat.2d.hdf5 ase_bedmachine_inverse_0_500m.2d.hdf5 post_inverse_0 post Cwshelf,muCoef,xVelb,yVelb,xVelo,yVelo,velc c_one,c_third,c_third_jreg_300, mu_coef,um,uo,uc


def post(c_one,mu_coef,u,v,uo,vo,uc,thk,topg,*etc):

    um = (u*u + v*v)**0.5
    uo = (uo*uo + vo*vo)**0.5
    c_third = c_one * (1+um)**(2.0/3.0)

    Nr = 100.0
    #ideally, min << 1 in the expression below but that
    #produces mad coeffcients at the present day GL.
    N = max(1, thk - max( -1027./917. * topg, 0) )
    
    
    #Joughin regularized rule
    uf = 300.0
    c_third_jreg_300 = c_third * ( um/uf + 1.0 )**(1.0/3.0) 

    #Joughin rule with pressure dependence
    #Tb ~ N |u|^1/3 / (Nr^3 * |u| + N^3 * |uf|)*1/3
    c_third_jreg_300_100 = c_third / N * (Nr**3* um/uf + N**3 )**(1.0/3.0)
    return c_one, c_third,  c_third_jreg_300, c_third_jreg_300_100,  mu_coef, um, uo, uc



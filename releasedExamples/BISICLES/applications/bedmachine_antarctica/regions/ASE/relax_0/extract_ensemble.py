#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:17:27 2021

@author: stephen
"""

import sys
import os
import glob
import numpy as np
sys.path.append(os.getcwd() + '/../../../python')
from bisiclesIO import BisiclesData, write_raster, read_raster, nctoamr
from ais_bedmachine import write_ais_nc

def opt_files():
    #list of the last file in any set of control runs that progressed beyonf iteration zero
    prefix = 'relax_muLT1_cont_mureg/ctrl.ase_bmach_muLT1_cont_mureg.03lev.'
    suffix = '.2d.hdf5'
    files = sorted(glob.glob('{}*[1-9]{}'.format(prefix,suffix)))
    opt_file = list() 
    ntime_prev = -1
    for f in files:
        ntime = f[len(prefix):len(prefix)+6]
        #niter = f[len(prefix)+6:len(prefix)+12]
        if (ntime == ntime_prev):
            opt_file[-1] = f
        else:
            opt_file.append(f)
        ntime_prev = ntime        
    return opt_file

def post(filename):
    #post processing for an hdf5 file
    #return the data defining an initial state
    print(filename)
    bike = BisiclesData(filename, level = 3, plot_file=False, iord=0)
   
    
    um = bike.speed

    c_third = bike.beta * (1.0 + um)**(2.0/3.0)
    c_third_jreg_300 = c_third * ( um/300.0 + 1.0 )**(1.0/3.0)
    c_third_jreg_600 = c_third * ( um/600.0 + 1.0 )**(1.0/3.0)
    c_third_jreg_1200 = c_third * ( um/1200.0 + 1.0 )**(1.0/3.0)
    
    mu_coef = bike.mucoef
    mu_coef_limit = (mu_coef - 1)*0.9 + 1
    
    #lower bedrock under floating ice for 'pure retreat' models.
    topg_retreat = np.where(bike.hab < 0.0, -1e+4, bike.topg)
    
    return bike.x,bike.y, {'thk':bike.thk, 'topg_retreat': topg_retreat, 'topg':bike.topg, 'mu_coef': mu_coef,
                           'mu_coef_limit': mu_coef_limit, 'misfit':bike.misfit,
            'c_one':bike.beta, 'c_third':c_third, 'c_third_jreg_300':c_third_jreg_300, 
            'c_third_jreg_600':c_third_jreg_600,'c_third_jreg_1200':c_third_jreg_1200}



files = opt_files()

for i,f in enumerate(files):
    if (i in [0,24]):
        x,y,var = post(f)   
        dx = x[1]-x[0]
        err = (np.sum((var['misfit'])**2))
        
        fn = 'u-ase-bmach-muLT1-ens-{:02d}.nc'.format(i)
        fh = 'u-ase-bmach-muLT1-ens-{:02d}.2d.hdf5'.format(i)
        write_ais_nc(fn,x,y,var)
        nctoamr(fn,fh, 'thk topg topg_retreat mu_coef mu_coef_limit c_one c_third c_third_jreg_300 c_third_jreg_600 c_third_jreg_1200' )

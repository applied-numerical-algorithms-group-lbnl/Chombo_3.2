#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 12:34:18 2020

@author: stephen
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

def read_pout(pout_file):
    
    csv_file = pout_file + '.csv'
    cmd = "awk '/grad/ {print $4 \",\" $14 \",\" $20}'  " + pout_file + " > " + csv_file
    print(cmd)
    os.system(cmd)
    df = pd.read_csv(csv_file,names=['fu','gc','gm'])
    return df
    

class L_curve_point:
    
    def __init__(self, misfit_norm, gradC_norm, gradMu_norm, name):
        self.misfit_norm = misfit_norm
        self.gradC_norm = gradC_norm
        self.gradMu_norm = gradMu_norm
        self.name = name


def read_pout_final(pout_file):
    
    p = read_pout(pout_file)
    n  = min(100,len(p))
        
    fu = p['fu'][n-1]
    gc = p['gc'][n-1]
    gm = p['gm'][n-1]
    
    return L_curve_point(fu,gc,gm, pout_file)
    



    
p104 = read_pout_final('regC1e4mu1e11/pout.ase_bmach.inverse..0')
p105 = read_pout_final('regC1e5mu1e12/pout.ase_bmach.inverse..0')
p102 = read_pout_final('regC1e2mu1e9/pout.ase_bmach.inverse..0')
p504 = read_pout_final('regC5e4mu5e11/pout.ase_bmach.inverse..0')
p304 = read_pout_final('regC3e4mu3e11/pout.ase_bmach.inverse..0')
p503 = read_pout_final('regC5e3mu5e10/pout.ase_bmach.inverse..0')
p153 = read_pout_final('regC15e3mu15e10/pout.ase_bmach.inverse..0')
p153_3010 = read_pout_final('regC15e3mu30e10/pout.ase_bmach.inverse..0')
p153_1010 = read_pout_final('regC15e3mu10e10/pout.ase_bmach.inverse..0')
p804 = read_pout_final('regC8e4mu8e11/pout.ase_bmach.inverse..0')
p204 = read_pout_final('regC2e4mu2e11/pout.ase_bmach.inverse..0')
for p in [p153,p104,p304,p504,p503,p804,p204,p153_3010]:
    plt.plot(p.misfit_norm,p.gradC_norm,'o',label=p.name[0:12])

plt.legend()
plt.ylabel(r'$|| \nabla C ||^2$')
plt.xlabel(r'$|| u - u_0 ||^2$')
#plt.xlim(2e14,5e14)

plt.savefig('lcurve_C.png',dpi=300)


plt.figure()
for p in [p153,p153_3010,p153_1010]:
    plt.plot(p.misfit_norm,p.gradMu_norm,'o',label=p.name)
plt.xlim(2e14,5e14)
plt.legend()
plt.ylabel(r'$|| \nabla \phi ||^2$')
plt.xlabel(r'$|| u - u_0 ||^2$')


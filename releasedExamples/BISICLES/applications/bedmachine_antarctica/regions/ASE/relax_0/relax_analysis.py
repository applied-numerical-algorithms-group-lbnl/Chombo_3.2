#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 10:42:55 2021

@author: stephen
"""

import numpy as np
import matplotlib.pyplot as plt
import re

cg = open('relax_with_no_sia/pout.ase_relax_0.0','r')

fm = list()

for cgline in cg:
    if re.match('CGOptimize iteration 0',cgline):
        fm.append(  float(cgline.split()[8]) )
        

plt.plot(fm,'b.-')
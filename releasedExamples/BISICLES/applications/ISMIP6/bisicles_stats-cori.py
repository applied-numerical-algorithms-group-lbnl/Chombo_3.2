#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:29:37 2015

@author: Dan

Call stats on plotfiles from a BISICLES run. Stores stats output in a plot*.hdf5.stats file for each plotfile

"""

import os
import csv
#import numpy as np
import glob
from netCDF4 import Dataset
from sys import argv
#from amrfile import io as amrio

script = argv[0]
stats = argv[1]
directory = argv[2]

def readstats(plotfile):
    
    #stats = "/scratch2/users/cashafer/MISMIP+/stats/stats2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex"
    statfile = str.replace(plotfile,'2d.hdf5','stats')
    nproc = '8'
    icedensity = '918.0'
    waterdensity = '1028.0'
    gravity = '9.81'
    
    if  ( not(os.path.isfile(statfile)) or 
            (os.path.getmtime(statfile) < os.path.getmtime(plotfile))): 
            
            cmd = "mpirun -np " + nproc + ' ' + stats + ' ' + plotfile + ' ' + icedensity + ' ' +  waterdensity + ' ' + gravity 
            print(cmd)
            os.system(cmd)
            #cmd2 = 'mv pout.0 ' + statfile
            cmd2 = 'grep time pout.0 > ' + statfile
            print(cmd2)
            os.system(cmd2)

    return statfile
            
def createstatsfile(name):
    
    print( 'scanning ' + name)
    pflist = sorted(glob.glob(name + "/plot.*.2d.hdf5"))

    outfile = name + '.allOut'
    cmd4 = 'touch ' + outfile
    os.system(cmd4)    
    cmd3 = 'mv -f ' + outfile  + ' save/' + outfile
    os.system(cmd3)
    os.system(cmd4)
    
    for plotfile in pflist:
        #main stats
        print (plotfile)
        statfile = readstats(plotfile)  
        print (statfile)

    cmd5 = 'cat ' + name + '/*.stats  > ' + outfile
    print (cmd5)
    os.system(cmd5)
       
    
names = glob.glob(directory)

for name in names:
    print (name)
    createstatsfile(name)

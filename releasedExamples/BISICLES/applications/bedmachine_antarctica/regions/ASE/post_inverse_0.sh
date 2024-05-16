flatten=$HOME/Development/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.OPT.ex
extract=$HOME/Development/BISICLES/code/filetools/extract2d.Linux.64.g++.gfortran.DEBUG.OPT.ex
pythonf=$HOME/Development/BISICLES/code/filetools/pythonf2d.Linux.64.g++.gfortran.DEBUG.OPT.ex

export PYTHONPATH=$PWD

$extract inverse_0/regC15e3mu15e10/ctrl.ase_bmach.03lev.000000000045.2d.hdf5 tmp.2d.hdf5 Cwshelf muCoef xVelb yVelb xVelo yVelo velc thickness Z_base
$flatten tmp.2d.hdf5 flat.2d.hdf5 3 # level 3: 500m
$pythonf flat.2d.hdf5 ase_bedmachine_inverse_0_500m.2d.hdf5 post_inverse_0 post Cwshelf,muCoef,xVelb,yVelb,xVelo,yVelo,velc,thickness,Z_base, c_one,c_third,c_third_jreg_300,c_third_jreg_300_100,mu_coef,um,uo,uc
rm tmp.2d.hdf5 flat.2d.hdf5

$extract inverse_0/regC15e3mu15e10_m10/ctrl.ase_bmach_muLT1.03lev.000000000063.2d.hdf5 tmp.2d.hdf5 Cwshelf muCoef xVelb yVelb xVelo yVelo velc
$flatten tmp.2d.hdf5 flat.2d.hdf5 3 # level 3: 500m
$pythonf flat.2d.hdf5 ase_bedmachine_inverse_0_muLT1_500m.2d.hdf5 post_inverse_0 post Cwshelf,muCoef,xVelb,yVelb,xVelo,yVelo,velc c_one,c_third,c_third_jreg_300,mu_coef,um,uo,uc
rm tmp.2d.hdf5 flat.2d.hdf5

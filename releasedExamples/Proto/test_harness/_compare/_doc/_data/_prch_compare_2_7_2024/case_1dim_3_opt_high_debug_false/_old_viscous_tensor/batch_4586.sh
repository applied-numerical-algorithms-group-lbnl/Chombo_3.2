#!/bin/csh

setenv CH_TIMER 1
setenv CH_OUTPUT_INTERVAL 1
setenv OMP_NUM_THREADS 1

setenv input case_1.inputs
setenv ProgName main.exe

module unload petsc hdf5
module load petsc/parallel hdf5/parallel

$ProgName $input >& pout.serial

#!/usr/bin/bash

module unload cmake
module unload cudatoolkit
module unload cray-hdf5-parallel
module unload craype-accel-nvidia80

module load cmake
module load cudatoolkit
module load cray-hdf5-parallel
module load craype-accel-nvidia80
alias allocstuff 'salloc --constraint=gpu --time=00:05:00 --cpus-per-task=10 --gpus=1 --qos=interactive --account=m1411_g'

#!/usr/bin/env bash 
# ======================================================================================
#             FILE: ismip6.sh
#
#            USAGE: ./ismip6.sh <name of experiment directory>
#
#      DESCRIPTION: Processes an experiment directory which is assumed to contain
#                   HDF5 plotfiles as well as HDF5 CF plotfiles contained within a CF
#                   directory. These plotfiles (both regular and CF) have variables that
#                   ismip6.sh processes individually. ismip6.sh then creates
#                   time-aggregated netCDF4 files for each variable of interest. 26 
#                   netCDF4 files are created for each of these 26 variables:
#		  
#                   lithk        topg        orog      base            xvelmean
#		    yvelmean     sftgif      acabf     libmassbffl     dlithkdt
#                   licalvf      ligroundf   sftgrf    sftflf          strbasemag
#                   lim          limnsw      iareag    iareaf          litemptop
#                   litempbotgr  litempbotfl tendacabf tendlibmassbffl tendlicalvf
#                   tendligroundf
#     
#         REQUIRES: - At least Python 2.7 (this code was tested with 2.7.6)
#		    - Up to date BISICLES version with working filetools
#			- glfacesNew2d
#			- extract2d
#			- flatten2d
#			- stats2d
#	    	    - Initial inputs file related to experiment being processed
#                   - 8km temperature field file
#                   - calculateTemperature.py
#		    - bisicles_stats-cori.py
#		    - calculateStatsfromAllOut.py
#		    - calculateFluxScalars.py
#                   - aggregate.py
#
# KNOWN WEAKNESSES: - Starting month/day/year is hard-coded in the following scripts and needs to
#                     be changed if you want your starting time to be correct in the
#                     final .nc files. The icesheet/institution/model name is also
#                     hard-coded in: 
#                       - calculateTemperature.py
#                       - calculateStatsfromAllOut.py
#                       - calculateFluxScalars.py
#                       - aggregate.py

#                   - Variables that aggregate.py looks to aggregate are hard-coded in
#                     lines 71, 106, 117, 129, and 143.
#
#                   - Expects that the CF files live within their own sub-directory
#                     called 'CF' inside the experiment directory that's being
#                     processed
#
#                   - When flattening, the value given depends on the total levels
#                     of the multi-level plotfile being flattened. In the code, level 0
#                     is assumed to be 8km (this value also appears in dependent 
#                     python scripts and may need to be changed if working with different
#                     leveled files) (we tested with level 3 plotfiles)
#
#                   - The code assumes that each plotfile in an experiment contains the
#                     expected variables and that none are missing. For example, if 
#                     plotfile #1, #2, #3 all have the sftgif, but plotfile #4 is
#                     missing, then the code will display a lot of errors, but will 
#                     continue to process.
#
#                   - Script can only handle square Antarctica grids
#
#                   - bisicles_stats-cori.py and calculateTemperature.py has ice
#                     density, water density, and gravity values that might differ from
#                     the values in the experiment that was ran.
#
#                   - This script assumes that the BISICLES executables are run using
#                     mpirun. If running in serial, search for "mpirun" and remove it.
#                     If running in parallel with something different, search for
#                     "mpirun" and edit.
#                      
#           AUTHOR: Courtney Shafer, cashafer201@gmail.com or cashafer@buffalo.edu
#
#             DATE: Sep 7th, 2021
# ======================================================================================

#           ============================
#             SET UP FOR USER (MODIFY)
#           ============================
# 0. Make sure the ISMIP6 processing directory is located somewhere that contains an
#    ample amount of space.
#
# 1. Create symbolic links to each of the experiment directories within the ismip6
#    processing directory
#	ex. Within the ismip6 processing directory (/ISMIP6), type 
#	    ln -s /location/of/an/experiment/directory/expt_1 expt_1
#	    to create the symbolic link "expt_1" that links to its original directory
# 
# 2. If your CF files are not already in their own CF directory within the experiment
#    directory, create the CF directory and move those over now.
#
# 3. Modify and uncomment variables below to match your own setup:
#
# Where is your BISICLES directory located (full directory name)?

BISICLES_HOME=/home/karloff/users/cashafer/bisicles/BISICLES

# Where is the ISMIP6 Processing directory located (full directory name)?

ISMIP6_HOME=/scratch2/users/cashafer/ISMIP6

# What are the names of the filetool executables? 
# Note: Make sure you are using the correct versions and change if needed
# For example: If you want MPI and Optimization enabled for your filetools, invoke:
# make glfacesNew2d OPT=TRUE MPI=TRUE
# when making the glfacesNew filetool (replace "glfacesNew" with the root name of
# the filetool that you're building)

GLFACESNEW=$BISICLES_HOME/code/filetools/glfacesNew2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex
EXTRACT=$BISICLES_HOME/code/filetools/extract2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex
FLATTEN=$BISICLES_HOME/code/filetools/flatten2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex
STATS=$BISICLES_HOME/code/filetools/stats2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex

# ismip6.sh needs the original inputs file used to generate the plotfiles. The inputs
# file is used to grab model parameters for glfacesNew. Copy it to the dependent_files
# directory then provide its name below.

INPUTS_FILE=$ISMIP6_HOME/dependent_files/inputs.example

# ismip6.sh needs an 8km resolution plotfile containing the temperature field of the
# ice sheet. The temperature field is used to generate 3 netCDF4 variables (litemptop,
# litempbotgr, litempbotfl) for submission. If the resolution of your temperature file
# is not in 8km, change the resolution of the file using the flatten filetool, and copy
# it to the dependent_files directory. Then, provide its name below.

TEMPERATURE_FIELD_FILE=$ISMIP6_HOME/dependent_files/example_antarctica-temperature-8km.2d.hdf5

# 4. Finally, to run, type ./ismip6.sh <name of linked experiment directory>
#                                     ** Do not put '/' in experiment name **
#
#           ==========================
#             SCRIPT (DO NOT MODIFY)
#           ==========================

# --------------------------------------------------------------------------------------
# I. DIRECTORY SETUP
# --------------------------------------------------------------------------------------

# Usage is ./ismip6.sh <linked experiment directory containing plotfiles>. If no
# experiment directory is provided, the code throws an error.

if [[ $# -eq 0 ]] ; then
    echo 'usage: ./ismip6.sh <name of experiment directory containing plotfiles>  ex: ./ismip6.sh ismip6_5 NOTE: Make sure the name of the experiment directory does not contain "/" at the end'
    exit 1
fi

EXPERIMENT=$1
echo Running $EXPERIMENT

# Within the ISMIP6 home directory /ISMIP6, directories necessary for processing are
# created for each experiment. Data that is created and/or extracted during processing
# will be saved in their appropriate experiment directory within these processing
# directories.

cd $ISMIP6_HOME

mkdir -p allOut_files/$EXPERIMENT
mkdir -p glfacesNew_output/$EXPERIMENT
mkdir -p extracted_variables/ligroundf/$EXPERIMENT
mkdir -p extracted_variables/sftgrf/$EXPERIMENT
mkdir -p extracted_variables/sftflf/$EXPERIMENT
mkdir -p extracted_variables/strbasemag/$EXPERIMENT
mkdir -p extracted_variables/thickness/$EXPERIMENT
mkdir -p extracted_variables/Z_base/$EXPERIMENT
mkdir -p extracted_variables/Z_bottom/$EXPERIMENT
mkdir -p extracted_variables/Z_surface/$EXPERIMENT
mkdir -p extracted_variables/xVel/$EXPERIMENT
mkdir -p extracted_variables/yVel/$EXPERIMENT
mkdir -p extracted_variables/iceFrac/$EXPERIMENT
mkdir -p extracted_variables/acabf/$EXPERIMENT/CF
mkdir -p extracted_variables/libmassbffl/$EXPERIMENT/CF
mkdir -p extracted_variables/dlithkdt/$EXPERIMENT/CF
mkdir -p extracted_variables/licalvf/$EXPERIMENT/CF
mkdir -p flattened_8km_ligroundf/$EXPERIMENT
mkdir -p flux_scalars/$EXPERIMENT
mkdir -p final_nc_files/$EXPERIMENT


# --------------------------------------------------------------------------------------
# II. EASY STATE VARIABLES ARE CALCULATED FIRST
# --------------------------------------------------------------------------------------

# Temperature data is first calculated using the calculateTemperature.py python script.
# calculateTemperature.py takes in the provided temperature field file and the initial
# plotfile of the experiment and calculates the surface temperature (litemptop),
# grounded basal temperature (litempbotgr), and the grounded floating temperature
# (litempbotfl). These are each calculated as single-level netCDF4 files and are saved
# in the appropriate experiment directory within the final_nc_files directory.  

echo Creating litemptop, litempbotgr, and litempbotfl .nc files from 8km temperature field file
cd $ISMIP6_HOME/final_nc_files/$EXPERIMENT

$ISMIP6_HOME/calculateTemperature.py $ISMIP6_HOME/$EXPERIMENT $TEMPERATURE_FIELD_FILE $EXPERIMENT

echo Done. litemptop, litempbotgr, and litempbotfl netCDF4 files located in final_nc_files/$EXPERIMENT

# Stats data is then calculated from running the bisicles_stats-cori.py and
# calculateStatsfromAllOut.py scripts. bisicles_stats-cori.py creates an .allOut file
# containing stats data for each timestep. calculateStatsfromAllOut.py goes through the
# .allOut file that was generated and grabs the data for the lim (icemassAll), limnsw
# (icemassAbove), iareag (groundedArea), and iareaf (floatingArea) variables and creates
# time-aggregated netCDF4 files for each variable which are saved in the appropriate
# experiment directory within final_nc_files.

echo Creating .stats files, .allOut file, and .nc files containing lim, limnsw, iareag, and iareaf

cd $ISMIP6_HOME

$ISMIP6_HOME/bisicles_stats-cori.py $STATS $EXPERIMENT
mv ${EXPERIMENT}.allOut $ISMIP6_HOME/allOut_files/$EXPERIMENT

cd $ISMIP6_HOME/final_nc_files/$EXPERIMENT

$ISMIP6_HOME/calculateStatsfromAllOut.py $ISMIP6_HOME $ISMIP6_HOME/allOut_files/$EXPERIMENT/${EXPERIMENT}.allOut $EXPERIMENT

echo Done. .allOut file is located in allOut_files/$EXPERIMENT and netCDF4 files are located in final_nc_files/$EXPERIMENT


# --------------------------------------------------------------------------------------
#  III. EXTRACTION OF VARIABLES FROM MULTI-COMPONENT FILES
# --------------------------------------------------------------------------------------

# Specific variables are contained within either the regular plotfiles or the CF files 
# and need to be extracted carefully. Some variables that are needed are not present in
# either filetypes and need to be derived. Specifically, glfacesNew takes in a regular 
# .hdf5 plotfile and derives the grounding line flux (ligroundf), grounded ice sheet
# area fraction (sftgrf), floating ice sheet area fraction (sftflf), and basal drag 
# (strbasemag). A new multi-level, multi-component .hdf5 file is created, which will
# then be used to extract the variables individually. This is performed first, then the
# remaining variables are extracted from the regular plotfiles and the CF files. All
# extracted variables are saved within the extracted_variables directory for further
# processing. 

# ==========================================================================
# |                    Variables extracted from files                      |
# ==========================================================================
# |____glfacesNew Output____|____regular plotfiles____|____CF plotfiles____|
# | ligroundf(F, FL)        | thickness(F, ST)        | acabf(F, FL)       |
# | sftgrf(F, ST)           | Z_base(F, ST)           | libmassbffl(F,FL)  |
# | sftflf(F, ST)           | Z_bottom(F, ST)         | dlithkdt(F, FL)    |
# | strbasemag(F, ST)       | Z_surface(F, ST)        | licalvf(F, FL)     |
# |                         | xVel(F, ST)             |                    |
# |                         | yVel(F, ST)             |                    |
# |                         | iceFrac(F, ST)          |                    |
# --------------------------------------------------------------------------


# glfacesNew output is generated first before extracting

cd $ISMIP6_HOME

echo Producing glfacesNew output

for i in $EXPERIMENT/plot*.hdf5
do
echo $i
mpirun -np 16 $GLFACESNEW $i $INPUTS_FILE glfacesNew_output/$i
done

echo Done - Data located in glfacesNew_output/$EXPERIMENT

# Then, each of the variables are extracted individually

cd $ISMIP6_HOME/glfacesNew_output

echo Extracting ligroundf, sftgrf, sftflf, strbasemag from glfacesNew output

for i in $EXPERIMENT/plot*.hdf5
do
echo $i
mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/ligroundf/$i ligroundf & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/sftgrf/$i sftgrf & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/sftflf/$i sftflf & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/strbasemag/$i strbasemag
done

echo Done extracting ligroundf, sftgrf, sftflf, and strbasemag. Data located in /extracted_variables/'<variable_name>'/$EXPERIMENT

# Next, the field state variables that live within the regular plotfiles are extracted
# individually 

cd $ISMIP6_HOME

echo Extracting thickness, Z_base, Z_surface, Z_bottom, xVel, yVel, and iceFrac from regular plotfiles

for i in $EXPERIMENT/plot*.hdf5
do
echo $i
mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/thickness/$i thickness & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/Z_base/$i Z_base & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/Z_bottom/$i Z_bottom & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/Z_surface/$i Z_surface & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/xVel/$i xVel & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/yVel/$i yVel & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/iceFrac/$i iceFrac
done
echo Done extracting thickness, Z_base, Z_surface, Z_bottom, xVel, yVel, and iceFrac. Data located in /extracted_variables/'<variable_name>'/$EXPERIMENT

# Finally the field flux variables that live within the CF files are extracted
# individually
cd $ISMIP6_HOME

echo Extracting acabf, libmassbffl, dlithkdt, and licalvf from CF plotfiles

for i in $EXPERIMENT/CF/plot*.hdf5
do
echo $i
mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/acabf/$i acabf & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/libmassbffl/$i libmassbffl & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/dlithkdt/$i dlithkdt & mpirun -np 16 $EXTRACT $i $ISMIP6_HOME/extracted_variables/licalvf/$i licalvf
done
echo Done extracting acabf, libmassbffl, dlithkdt, and licalvf. Data located in /extracted_variables/'<variable_name>'/$EXPERIMENT

# --------------------------------------------------------------------------------------
# IV. SPATIAL INTEGRATION TO PRODUCE FLUX SCALAR VARIABLES
# --------------------------------------------------------------------------------------

# There are a few remaining variables that need to be derived, specifically the flux
# scalar variables, which are total SMB flux (tendacabf), total BMB flux beneath
# floating ice (tendlibmassbffl), total calving flux (tendlicalvf), and total grounding
# line flux (tendligroundf). calculateFluxScalars.py takes in a single-level file and
# performs the spatial integration to produce these variables. The CF files themselves
# which contain the first three variables are already single-level, however, the
# ligroundf files that were extracted earlier are not, so they need to be flattened
# first.     

# First, the extracted ligroundf files are flattened to level 0 (8km resolution)
cd $ISMIP6_HOME/extracted_variables/ligroundf

echo Flattening extracted ligroundf plotfiles 

for i in $EXPERIMENT/plot*.hdf5
do
echo $i
mpirun -np 16 $FLATTEN $i $ISMIP6_HOME/flattened_8km_ligroundf/$i 0
done
echo Done flattening ligroundf. Data located in /flattened_8km_ligroundf/$EXPERIMENT

# Then the flux scalars are calculated. calculateFluxScalars.py takes in the ligroundf
# directory location as well as the CF files directory location to calculate these
# simultaneously.

cd $ISMIP6_HOME/flux_scalars/$EXPERIMENT

echo Calculating the flux scalars
$ISMIP6_HOME/calculateFluxScalars.py $ISMIP6_HOME/flattened_8km_ligroundf/$EXPERIMENT $ISMIP6_HOME/$EXPERIMENT/CF $EXPERIMENT
echo Done calculating the flux scalars. Data located in /flux_scalars/$EXPERIMENT

cp -r $ISMIP6_HOME/flux_scalars/$EXPERIMENT/. $ISMIP6_HOME/final_nc_files/$EXPERIMENT

# --------------------------------------------------------------------------------------
# V. AGGREGATING REMAINING HDF5 FILES FOR EACH VARIABLE INTO SINGLE NETCDF4 FILE
# --------------------------------------------------------------------------------------

# The variables that were extracted in section III still need to converted into a single
# netCDF4 file. aggregate.py takes in these files for each variable and aggregates them
# into a single netCDF4 file. It can handle both single-level and multi-level data, but
# it is assumed that level 0 is 8km. The variables are hard-coded in the code. 

echo Aggregating plotfiles into final netCDF4 files
cd $ISMIP6_HOME/final_nc_files/$EXPERIMENT

$ISMIP6_HOME/aggregate.py $EXPERIMENT $ISMIP6_HOME/extracted_variables
echo Done.

 


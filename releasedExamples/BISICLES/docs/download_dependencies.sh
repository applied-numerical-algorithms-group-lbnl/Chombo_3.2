cd $BISICLES_HOME
echo `pwd`



#get hdf5 sources
if !(test -e hdf5-1.10.9.tar.bz2) then
    echo "downloading hdf5"
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.9/src/hdf5-1.10.9.tar.bz2
fi
mkdir -p hdf5/parallel/src
tar -jxf  hdf5-1.10.9.tar.bz2 -C hdf5/parallel/src

mkdir -p hdf5/serial/src
tar -jxf  hdf5-1.10.9.tar.bz2 -C hdf5/serial/src


#get netcdf sources

if !(test -e netcdf-4.1.3.tar.gz) then
    echo "downloading netcdf-4.1.3.tar.gz"
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.1.3.tar.gz
fi
mkdir -p netcdf/parallel/src
tar -zxf netcdf-4.1.3.tar.gz -C netcdf/parallel/src

mkdir -p netcdf/serial/src
tar -zxf netcdf-4.1.3.tar.gz -C netcdf/serial/src





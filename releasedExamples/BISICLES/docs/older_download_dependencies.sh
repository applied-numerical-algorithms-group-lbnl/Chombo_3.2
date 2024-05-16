cd $BISICLES_HOME
echo `pwd`

#get hdf5 sources
if !(test -e hdf5-1.6.10.tar.gz) then
    echo "downloading hdf5-1.6.10.tar.gz"
    wget http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.6.10/src/hdf5-1.6.10.tar.gz
fi
mkdir -p hdf5/parallel/src
tar -zxf hdf5-1.6.10.tar.gz -C hdf5/parallel/src

mkdir -p hdf5/serial/src
tar -zxf hdf5-1.6.10.tar.gz -C hdf5/serial/src


#get hdf5 sources
#if !(test -e hdf5-1.8.9.tar.bz2) then
#    echo "downloading hdf5-1.8.9.tar.gz"
#    wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.9.tar.bz2
#fi
#mkdir -p hdf5/parallel/src
#tar -jxf  hdf5-1.8.9.tar.bz2 -C hdf5/parallel/src

#mkdir -p hdf5/serial/src
#tar -jxf  hdf5-1.8.9.tar.bz2 -C hdf5/serial/src


#get netcdf sources

if !(test -e netcdf-4.1.2.tar.gz) then
    echo "downloading netcdf-4.1.2.tar.gz"
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.1.2.tar.gz
fi
mkdir -p netcdf/parallel/src
tar -zxf netcdf-4.1.2.tar.gz -C netcdf/parallel/src

mkdir -p netcdf/serial/src
tar -zxf netcdf-4.1.2.tar.gz -C netcdf/serial/src





# Chombo_3_2/example/TimeParallel

## Tools for time parallelism using MPI communicator splitting.
* BoxIterator owns a communicator.
* Most of Chombo does not care.
* HDF5 and BoxLayoutData do care.

## Directories
* comm_split mimics an MPI tutorial to split grids into multiple colors.
* split_file_IO shows how do input and output in the context of colors.
* coloredSolves shows how the whole dance would work.

## Makefile stuff
* None of this makes any sense if MPI=FALSE.
* HDF5 needs to be configured for parallel. 


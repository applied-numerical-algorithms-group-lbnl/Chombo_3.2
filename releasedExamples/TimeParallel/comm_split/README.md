# comm_split_test
* Parallel in time algorithms require MPI communicators to be split.
* The solution  data is replicated a specfied number of times.
* We call MPI_COMM_SPLIT to get a number of communicators
* The way this works:
  * MPI_COMM_WORLD is the same as Chombo_MPI::comm
  * We say we want to have N replications of the solution data.
  * MPI_Comm_Split is used to make N communicators (called row communicators).
  * Each processor will know about the world communcator and one row comunicator.
  * Each process knows its place in both the world and its row communcators.
  * Data transfer for this experimental bit will be via file I/O using HDF5.

# comm_split_test
* Parallel in time algorithms require MPI communicators to be split.
* The solution  data is replicated a specfied number of times.
* We call MPI_COMM_SPLIT to get a number of communicators
* This particular test runs AMRPoisson in three
  * First it runs the standard chombo way to make sure all that still works.
  * Second, you pick a number of colors (N) and it runs the solver N times.
  * Third it runs the solver N times simultaneously.
    *  The communicator is split into N colors.
    *  Each proc sees the global solve and its color solve.
  * All this stuff gets timed.
  * This test only makes sense if compiled with MPI=TRUE so there may be missing guards.

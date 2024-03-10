# Chombo_3.2/releasedExamples/Proto/test_harness/_mpi_strong
* There are a bunch of examples and a bunch of input files and both 2 and 3 dimensions.
* For this test,  there are also varying numbers of processors.
* This test runs all these configurations puts results in managable places.
* The lofty purpose here was just to make sure everything worked in parallel.
* We also get some sense of the strong scaling behavior of each operator.
* As of 3-10-2024, everything runs in parallel.
* See document in ../../documents/mpi_strong.pdf.
* Many directories are created.
* A single script to run them all is also created.


# Script configure_strong_scaling.py: 
* makes many directories
* compiles everything
* puts batch scripts and input files in the  right places.
* makes one script to rule them all and in darkness bind them (dirname/runall.sh).

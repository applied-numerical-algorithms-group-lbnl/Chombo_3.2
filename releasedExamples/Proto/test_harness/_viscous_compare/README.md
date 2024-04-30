# Chombo_3.2/releasedExamples/Proto/test_harness/_viscous_compare
* There are a 3 different viscous operators and a bunch of input files and both 2 and 3 dimensions.
* For this test,  there are also varying numbers of processors.
* This test runs all these configurations puts results in managable places.
* The purpose here is to compare the old, fortran-based operators with two different viscous operators.
* We also get some sense of the strong scaling behavior of each operator.
* As of 4-30-2024, everything runs in parallel.
* Many directories are created.
* A single script to run them all is also created.


# Script configure_viscous_compare.py
* makes many directories
* compiles everything
* puts batch scripts and input files in the  right places.
* makes one script to rule them all and in darkness bind them (dirname/runall.sh).

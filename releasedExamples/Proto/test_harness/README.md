# Chombo_3.2/releasedExamples/Proto/test_harness
* Some simple python-driven utilities to set up and run the PrCh examples.
* Simple macros are put into an input template.
* You can set up and compile a bunch of differrent configurations (DIM, MPI, etc) at once.
* Many directories are created.
* A single script to run them all is also created.
* For supercomputers, this means a bunch of batch submissions.

# Directories
* _input_templates holds templated versions of the input files for the Proto/Chombo examples.
* _batch_templates are batch scripts to run on the platforms I use.
* _compare runs the examples in serial and puts the results in manageable places.   See ../documents/compare.pdf.
* _mpi_strong tests the strong scaling  behavior of each example and puts the results in managable place.  It includes a script for data retrieval.  See ../documents/mpi_strong.pdf.
* _viscous_compare compares the three diffeent viscous operators (fortran-based and two proto discretizations).


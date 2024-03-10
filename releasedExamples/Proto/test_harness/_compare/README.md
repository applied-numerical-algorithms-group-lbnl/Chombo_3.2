# Chombo_3.2/releasedExamples/Proto/test_harness/compare
* There are 8 examples (4 old, 4 new), 4 cases (input files) and both two and three dimensions. 
* This test runs all 64 of these configurations (in serial) and puts results in managable places.
* The lofty purpose here was just to make sure everything worked.
* Using these comparisons, I was able to fix some multigrid performance issues.
* As of 3-10-2024 both old and new  perform similarly.
* See document in ../../documents/compare.pdf.
* Many directories are created.
* A single script to run them all is also created.


# Contents
* _doc documentation source files
* configure_test Python script that:
** makes many directories
** compiles everything
** puts batch scripts and input files in the  right places.
** makes one script to rule them all and in darkness bind them (dirname/runall.sh).

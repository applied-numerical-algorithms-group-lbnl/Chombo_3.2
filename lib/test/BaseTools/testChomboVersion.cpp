#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <string>

#include "CHOMBO_VERSION.H"
#include "SPMD.H"
#include "parstream.H"

#include "UsingBaseNamespace.H"

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  string testVersion;
  string testFlags;
  string testDate;
  string testTime;

  getChomboVersion(testVersion,testFlags,testDate,testTime);

  pout() << "Version: " << testVersion << endl;
  pout() << "Flags:   " << testFlags   << endl;
  pout() << "Date:    " << testDate    << endl;
  pout() << "Time:    " << testTime    << endl;

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}


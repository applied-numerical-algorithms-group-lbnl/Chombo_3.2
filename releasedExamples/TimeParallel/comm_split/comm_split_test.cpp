#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "memusage.H"

int spmdTest()
{
  return 0;
}

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;

  // scoping...
  {
    if (argc < 2)
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int solverStatus = spmdTest();
    status += solverStatus;
  }
  //end scoping trick
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}

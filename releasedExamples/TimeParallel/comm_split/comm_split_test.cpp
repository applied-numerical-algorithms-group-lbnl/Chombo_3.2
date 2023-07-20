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

/**
   just to get something working
   taken unedited from mpitutorial.com
**/
int
mpiTutorialCommSplitTest()
{
  // Get the rank and size in the original communicator
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int color = world_rank / 4; // Determine color based on row

// Split the communicator based on the color and use the
// original rank for ordering
  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

  int row_rank, row_size;
  MPI_Comm_rank(row_comm, &row_rank);
  MPI_Comm_size(row_comm, &row_size);

//  printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",
//         world_rank, world_size, row_rank, row_size);
  //chombofied output so this can be run at large concurrencies.
  pout() << "mpi_tutorial: world rank = " << world_rank << endl;
  pout() << "mpi_tutorial: world size = " << world_size << endl;
  pout() << "mpi_tutorial: row   rank = " <<   row_rank << endl;
  pout() << "mpi_tutorial: row   size = " <<   row_size << endl;
  MPI_Comm_free(&row_comm);

  return 0;
}
///
/**
   Taking the above tutorial and chombofying it.
 **/
int
chomboCommSplitTest()
{
  ParmParse pp("chomboCommSplitTest");
  int nx, maxBoxSize, numSubset;
  pp.get("nx", nx);
  pp.get("maxBoxSize", maxBoxSize);
  pp.get("numSubset" , numSubset);
  ///
  // Get the rank and size in the original communicator
  int world_rank, world_size;
  MPI_Comm_rank(Chombo_MPI::comm, &world_rank);
  MPI_Comm_size(Chombo_MPI::comm, &world_size);

  // Determine color based on original rank
  int icolor = world_rank / numSubset; 

  MPI_Comm row_comm;
  MPI_Comm_split(Chombo_MPI::comm, icolor, world_rank, &row_comm);
  
  int row_rank, row_size;
  MPI_Comm_rank(row_comm, &row_rank);
  MPI_Comm_size(row_comm, &row_size);
  
  pout() << "chomboCommSplit: world rank = " << world_rank << endl;
  pout() << "chomboCommSplit: world size = " << world_size << endl;
  pout() << "chomboCommSplit: row   rank = " <<   row_rank << endl;
  pout() << "chomboCommSplit: row   size = " <<   row_size << endl;
  pout() << "chomboCommSplit: proc color = " << icolor     << endl;
  
  IntVect ivlo =        IntVect::Zero;
  IntVect ivhi = (nx-1)*IntVect::Unit;
  Box domain(ivlo, ivhi);
  Vector<Box> boxes;
  Vector<int> procs;
  domainSplit(domain, boxes,  maxBoxSize);
  LoadBalance(procs, boxes);
  DisjointBoxLayout(boxes, procs);
    
  
  return 0;
}

/// drives the tests here
int spmdTest()
{
  {
    int ierr = mpiTutorialCommSplitTest();
    if( ierr != 0)
    {
      cerr << " mpiTutorialCommSplitTest retuned error code "  << ierr;
      return ierr;
    }
  }
  { 
    int ierr = chomboCommSplitTest();
    if( ierr != 0)
    {
      cerr << " mpiTutorialCommSplitTest retuned error code "  << ierr;
      return ierr;
    }
  }
  return 0;
}

/// init mpi, init parmparse.  call test.  dump timers. finalize mpi.
int main(int argc, char* argv[])
{
  ///init mpi
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  if (argc < 2)
  {
    cerr << " usage " << argv[0] << " <input_file_name> " << endl;
    exit(0);
  }

  ///init parmparse using the ancient formula
  char* in_file = argv[1];
  ParmParse  pp(argc-2,argv+2,NULL,in_file);

  /// call test
  int status = spmdTest();
  if(status != 0)
  {
    cerr << "spmdTest returned error code " << status << endl;
    return status;
  }
  
#ifdef CH_MPI
  dumpmemoryatexit();
  ///dump timers
  CH_TIMER_REPORT();
  ///finalize mpi
  MPI_Finalize();
#endif

  return(0);
}

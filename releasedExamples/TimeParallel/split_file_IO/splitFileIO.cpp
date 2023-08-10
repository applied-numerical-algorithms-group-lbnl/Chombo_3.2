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
#include "BCFunc.H"
#include "AMRPoissonOp.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CH_HDF5.H"




int
writeAndReadData()
{
  ParmParse pp("writeAndReadData");
  int nx, maxBoxSize, procsPerColor;
  pp.get("nx"        , nx         );
  pp.get("maxBoxSize", maxBoxSize );
  pp.get("procsPerColor" , procsPerColor  );

  pout() << "writeAndReadData:nx         =" << nx         << endl;
  pout() << "writeAndReadData:maxBoxSize =" << maxBoxSize << endl;
  pout() << "writeAndReadData:procsPerColor  =" << procsPerColor  << endl;
  IntVect ivlo =        IntVect::Zero;
  IntVect ivhi = (nx-1)*IntVect::Unit;
  Box domain(ivlo, ivhi);
  Vector<Box> boxes;
  Vector<int> procs;
  domainSplit(domain, boxes,  maxBoxSize);
  double dx = 1./nx;
  LoadBalance(procs, boxes);

  ///
  // Get the rank and size in the original communicator
  int world_rank, world_size;
  MPI_Comm_rank(Chombo_MPI::comm, &world_rank);
  MPI_Comm_size(Chombo_MPI::comm, &world_size);

  CH_TIME("writeAndReadData: actual test");
  // Determine color based on original rank
  int icolor = world_rank / procsPerColor;

  MPI_Comm color_comm;
  MPI_Comm_split(Chombo_MPI::comm, icolor, world_rank, &color_comm);
  
  int color_rank, color_size;
  MPI_Comm_rank(color_comm, &color_rank);
  MPI_Comm_size(color_comm, &color_size);
  
  pout() << "writeAndReadData: world rank = " << world_rank << endl;
  pout() << "writeAndReadData: world size = " << world_size << endl;
  pout() << "writeAndReadData: color rank = " << color_rank << endl;
  pout() << "writeAndReadData: color size = " << color_size << endl;
  pout() << "writeAndReadData: proc color = " << icolor     << endl;
  

  DisjointBoxLayout dblWorld(boxes, procs, Chombo_MPI::comm);
  DisjointBoxLayout dblColor(boxes, procs, color_comm);
  writeColorFileFromColor(dblColor, domain, dx,   icolor);
  readColorFileIntoWorld( dblWorld, domain, dx,   icolor);

  return 0;
}


/// init mpi, init parmparse.  call test.  dump timers. finalize mpi.
int main(int argc, char* argv[])
{
  ///init mpi
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 4586;
  if (argc < 2)
  {
    cerr << " usage " << argv[0] << " <input_file_name> " << endl;
    return status;
  }

  ///init parmparse using the ancient formula
  char* in_file = argv[1];
  ParmParse  pp(argc-2,argv+2,NULL,in_file);

  /// call test
  status = writeAndReadData();
  if(status != 0)
  {
    cerr << "writeAndReadData returned error code " << status << endl;
    return status;
  }
  
#ifdef CH_MPI
  dumpmemoryatexit();
  ///dump timers
  CH_TIMER_REPORT();
  ///finalize mpi
  MPI_Finalize();
#endif

  return(status);
}

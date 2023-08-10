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

////
int
writeColorFileFromColor(const DisjointBoxLayout  & a_grid,
                        const ProblemDomain      & a_domain,
                        const Real               & a_dx,
                        int       a_icolor) //included in file names if positive
{
  string phi_fname("4586.hdf5");
  //negative signals world
  if(a_icolor >= 0)
  {
    phi_fname = string("phi_color") + to_string(a_icolor) + string("_chombo.hdf5");
  }
  else
  {
    pout() << "writeColorFromColor: invalid color" << endl;
    return -4586;
  }
  ///packaging up stuff so I can call the most general write function
  Vector<DisjointBoxLayout> amrGrids(1, a_grid);
  Vector<ProblemDomain>   amrDomains(1, a_domain);
  Vector<Real> amrDx(1, a_dx);
  Vector<int> refRatios(1, 2);

  double time = 0;  double dt = 0; //just used for labeling
  Vector<string> phi_vars(1, string("phi"));
  Box dombox= a_domain.domainBox();
  LevelData<FArrayBox> philev(a_grid, 1, IntVect::Zero);
  DataIterator dit = a_grid.dataIterator();
  //have to set it to something
  Real value = 4586/Real(a_icolor + 1);
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    philev[dit[ibox]].setVal(value);
  }
  Vector<LevelData<FArrayBox> * > phi(1, &philev);
  pout() << "writeColorFromColor:  writing phi to " << phi_fname << endl;
  WriteAMRHierarchyHDF5(phi_fname, amrGrids, phi, phi_vars, dombox, a_dx, dt, time, refRatios, 1);

  return 0;
}

////
int
readColorFileIntoWorld(int       a_icolor) 
{
  string phi_fname("4586.hdf5");
  //negative signals world
  if(a_icolor >= 0)
  {
    phi_fname = string("phi_color") + to_string(a_icolor) + string("_chombo.hdf5");
  }
  else
  {
    pout() << "readColorIntoWorld: invalid color" << endl;
    return -4586;
  }
  pout() << "readColorIntoWorld: reading" << endl;

  Vector<DisjointBoxLayout>       grids;
  Vector<LevelData<FArrayBox>* >  phi;
  Vector<string>                  var_names;
  Box                             domain;
  Real                            dx, dt, time;
  Vector<int>                     refRatio;
  int                             numLevels;

  ReadAMRHierarchyHDF5( phi_fname,
                        grids,    phi, var_names, domain,   
                        dx, dt, time, refRatio,  numLevels);
                     
  return 0;
}


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
  MPI_Comm_size(color_comm, &color_size);              // 
  
  pout() << "writeAndReadData: world rank = " << world_rank << endl;
  pout() << "writeAndReadData: world size = " << world_size << endl;
  pout() << "writeAndReadData: color rank = " << color_rank << endl;
  pout() << "writeAndReadData: color size = " << color_size << endl;
  pout() << "writeAndReadData: proc color = " << icolor     << endl;
  

  DisjointBoxLayout dblWorld(boxes, procs, Chombo_MPI::comm);
  DisjointBoxLayout dblColor(boxes, procs, color_comm);
  {
    int eekwrite = writeColorFileFromColor(dblColor, domain, dx,   icolor);
    if(!eekwrite) return eekwrite;
  }
  {
    //read makes its own dbl
    int eekread  = readColorFileIntoWorld(icolor);
    if(!eekread) return eekread;
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

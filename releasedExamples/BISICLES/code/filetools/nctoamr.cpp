 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// nctoamr.cpp
// read in data from a netcdf file, write out a single level AMR hierarchy
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "fabncio.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"


int main(int argc, char* argv[]) {

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif
  

    if(argc < 4) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> <var 1> [<var 2>, ...] " 
		  << std::endl; 
	exit(0); 
      }
    
    char* in_file = argv[1];
    char* out_file = argv[2];
    
    Vector<std::string> var;
    for (int i = 3; i < argc; i++)
      {
	var.push_back(std::string(argv[i]));
      }  
    
    //pout() << "converting netcdf file " << in_file << " to AMR file " << out_file << std::endl;
    
    Box box;
    Real dx;
    DisjointBoxLayout grids;
    LevelData<FArrayBox> levelData;
    
    // scope for copying into a distributed holder
    { 
      FArrayBox fab;
      if (procID() == uniqueProc(SerialTask::compute))
        {
#ifdef HAVE_NETCDF
          NCIO::readFAB(in_file,var,fab,dx);
	  box = fab.box();
#else
          MayDay::Error("netcdf input requested but netcdf support not built");
#endif
	    }// end if serial compute
      
      broadcast(dx, uniqueProc(SerialTask::compute));
      broadcast(box,uniqueProc(SerialTask::compute));
      ProblemDomain pd(box);
      Vector<Box> boxes(1,pd.domainBox());
      Vector<int> procAssign(1,uniqueProc(SerialTask::compute));
      DisjointBoxLayout serialGrids(boxes, procAssign, pd);
      LevelData<FArrayBox> serialLevelData(serialGrids,var.size(),IntVect::Zero);
      
      for (DataIterator dit(serialGrids); dit.ok(); ++dit)
	{
	  serialLevelData[dit].copy(fab,0,0,fab.nComp());
	}
      
      // now distribute if needed
      // arbitrarily pick 1024 as max box size
      int maxBoxSize = 1024;
      domainSplit(box, boxes, maxBoxSize);
      LoadBalance(procAssign, boxes);
      grids.define(boxes, procAssign, pd);
      levelData.define(grids, var.size(), IntVect::Zero);
      serialLevelData.copyTo(levelData);
      
    }
    
    Vector<LevelData<FArrayBox>* > vectData(1,&levelData);
    Vector<DisjointBoxLayout > vectGrids(1,grids);
    Vector<int > vectRatio;
    
    const Real dt = 1.0;
    const Real time = 0.0;
    WriteAMRHierarchyHDF5(std::string(out_file), vectGrids, vectData, var , 
			  box, dx, dt, time, vectRatio, 1);
    
    
  
    
  }  // end nested scope
  CH_TIMER_REPORT();
  
#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

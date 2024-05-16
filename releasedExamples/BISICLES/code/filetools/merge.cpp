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
// merge.cpp
// read in two multi-level AMR hierarchies and write their union.
// prefix component names from the first "a." and the second "b."
//===========================================================================
#include <iostream>
#include "AMRIO.H"
#include "NamespaceHeader.H"
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

    if(argc != 4) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file a> <input_file b> <output_file>" << std::endl; 
	exit(0); 
      }

    char* in_file_a = argv[1];
    char* in_file_b = argv[2];
    char* out_file = argv[3];


    Vector<std::string> names_a;
    Vector<LevelData<FArrayBox>* > data_a;
    Vector<DisjointBoxLayout> grids_a;
    Vector<int> ratio_a;
    Real crseDx_a = 0.0, dt_a = 0.0, time_a = 0.0;
    Box crseBox_a;
    int numLevels_a;
    int status = ReadAMRHierarchyHDF5
      (in_file_a,grids_a,data_a,names_a,crseBox_a,crseDx_a,dt_a,time_a,
       ratio_a,numLevels_a);
    if (status != 0)
      {
	MayDay::Error("failed to read AMR hierarchy a");
      }
    
    Vector<std::string> names_b;
    Vector<LevelData<FArrayBox>* > data_b;
    Vector<DisjointBoxLayout> grids_b;
    Vector<int> ratio_b;
    Real crseDx_b = 0.0, dt_b = 0.0, time_b = 0.0;
    Box crseBox_b;
    int numLevels_b;
    status = ReadAMRHierarchyHDF5
      (in_file_b,grids_b,data_b,names_b,crseBox_b,crseDx_b,dt_b,time_b,
       ratio_b,numLevels_b);
    if (status != 0)
      {
	MayDay::Error("failed to read AMR hierarchy b");
      }

    if (numLevels_b != numLevels_a)
      {
	MayDay::Error("numLevels_b != numLevels_a");
      }

    if (crseDx_b != crseDx_a)
      {
	MayDay::Error("crseDx_b != crseDx_a");
      }

    Vector<std::string> names_merge(names_a.size() + names_b.size());
    for (int i = 0; i < names_a.size(); i++)
      {
	names_merge[i] = names_a[i] + ".a";
      }
    for (int i = 0; i < names_b.size(); i++)
      {
	names_merge[i + names_a.size()] = names_b[i] + ".b";
      }

    Vector<LevelData<FArrayBox>* > outData(numLevels_a,NULL);
    for (int lev = 0; lev < numLevels_a; lev++)
      {
	if (ratio_b[lev] != ratio_a[lev])
	  {
	    MayDay::Error("ratio_b[lev] != ratio_a[lev]");
	  }
	outData[lev] = new LevelData<FArrayBox>(grids_a[lev],names_merge.size(),data_a[lev]->ghostVect());
	DataIterator dit_a(grids_a[lev]);
	DataIterator dit_b(grids_b[lev]);
	for (dit_a.reset(),dit_b.reset();dit_a.ok(),dit_b.ok();++dit_a,++dit_b)
	  {
	    (*outData[lev])[dit_a].copy( (*data_a[lev])[dit_a], 0,0,names_a.size());
	    (*outData[lev])[dit_a].copy( (*data_b[lev])[dit_b], 0,names_a.size(),names_b.size());
	  }
      }

    WriteAMRHierarchyHDF5(out_file, grids_a, outData, names_merge, 
			  crseBox_a, crseDx_a, dt_a, time_a, ratio_a,numLevels_a);
    
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}

#include "NamespaceFooter.H"

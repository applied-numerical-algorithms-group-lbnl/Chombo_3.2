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
// rescale.cpp
//===========================================================================

#include <iostream>
#include "AMRIO.H"
#include "fabncio.H"
#include "NamespaceHeader.H"

bool verbose = true;

enum out_file_type_enum {hdf5,nc};

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
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> scale" << std::endl; 
	exit(0); 
      }

    char* in_file = argv[1];
    char* out_file = argv[2];

 
    out_file_type_enum out_file_type;
    out_file_type = hdf5;

    Real scale = atof(argv[3]);

    if (verbose)
      {
        pout ();
        pout() << "reading AMR file..." << endl;
      }
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout> grids;
    Vector<string> names;
    Vector<int> ratio;
    Real crseDx = 0.0, dt = 0.0, time = 0.0;
    Box crseBox;
    int numLevels;
    int status = ReadAMRHierarchyHDF5
      (in_file,grids,data,names,crseBox,crseDx,dt,time,
       ratio,numLevels);
    if (status != 0)
      {
        MayDay::Error("failed to read AMR hierarchy");
      }
    
    if (verbose)
      {
        pout() << "... done." << endl;
      }

    for (int lev=0; lev<data.size(); lev++)
      {
        LevelData<FArrayBox>& dataLev = *data[lev];
        DataIterator dit = dataLev.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            dataLev[dit].mult(scale);
          }
      }

    WriteAMRHierarchyHDF5(out_file, grids, data, names, 
                          crseBox, crseDx, dt, time, ratio, numLevels);
    
    // clean up memory
    for (int lev=0; lev<data.size(); lev++)
      {
        if (data[lev] != NULL)
          {
            delete data[lev];
            data[lev] = NULL;
          }
      }
    
  }  // end nested scope
  CH_TIMER_REPORT();
  
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}
#include "NamespaceFooter.H"

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
// extract.cpp
// read a multi-level AMR hierarchy, select a subset of its components, and
// write the result to a Chombo hdf5 file
//===========================================================================
#include <iostream>
#include "AMRIO.H"
#include "NamespaceHeader.H"

bool verbose = true;

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
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> <var 1> [<var 2> [...]]" << std::endl; 
	exit(0); 
      }

    char* in_file = argv[1];
    char* out_file = argv[2];

    Vector<std::string> var;
    for (int i = 3; i < argc; i++)
      {
	var.push_back(std::string(argv[i]));
      }  

    if (verbose) 
      {
        pout() << "reading " << in_file << "..." << endl;
      }

    Vector<std::string> names;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout> grids;
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
        pout() << endl;
        pout() << "looking for variables to extract..." << endl;
      }
    
    
    Vector<int> outToIn(var.size());
      
    for (int i =0; i < var.size(); i++)
      {
	bool found = false;
	for (int j =0; j < names.size(); j++)
	  {
            found = (names[j] == var[i]);
	    if (found)
	      {
		outToIn[i] = j;
		break;
	      }
	  }
	if (!found)
	  {
	    std::cerr << var[i] << " not found " << std::endl;
	    MayDay::Error("variable not found");
	  }
      }
    
    if (verbose)
      {
        pout() << "... done." << endl;
        pout() << "constructing output data" << endl;
      }

    Vector<LevelData<FArrayBox>* > outData(numLevels,NULL);
    for (int lev = 0; lev < numLevels; lev++)
      {
	
	outData[lev] = new LevelData<FArrayBox>(grids[lev],var.size(),data[lev]->ghostVect());
	for (int i = 0 ;i < outToIn.size(); i++)
	  {
	    //data[lev]->copyTo(Interval(outToIn[i],outToIn[i]),*outData[lev],Interval(i,i));
	    for (DataIterator dit(grids[lev]);dit.ok();++dit)
	      {
		(*outData[lev])[dit].copy( (*data[lev])[dit], outToIn[i], i, 1);
	      }
	  }
      }

    if (verbose)
      {
        pout() << "writing file..." << endl;
      }
    WriteAMRHierarchyHDF5(out_file, grids, outData, var, 
			  crseBox, crseDx, dt, time, ratio,numLevels);
    

    if (verbose)
      {
        pout() << "... done." << endl;
      }
  
    for (int lev=0; lev<data.size(); lev++)
      {
        if (data[lev] != NULL)
          {
            delete data[lev];
            data[lev] = NULL;
          }
        if (outData[lev] != NULL)
          {
            delete outData[lev];
            outData[lev] = NULL;
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
